import type { ShikiTransformer, ThemedToken } from 'shiki'
import type { Element } from 'hast'
import fs from 'node:fs'
import path from 'node:path'

/*
 * Links identifiers in Julia code blocks to their docstrings and attaches the data
 * for a hover tooltip. This is the VitePress stand-in for the reference links and
 * popups of DocumenterCodeBlocks.jl, which only post-processes Documenter.HTML output.
 *
 * The symbol index is read from the markdown DocumenterVitepress emits into
 * `build/.documenter`. That markdown is complete before `vitepress build` starts,
 * so this works both when the build runs inside `makedocs` (CI) and when make.jl
 * runs it afterwards (Windows).
 *
 * Unlike DocumenterCodeBlocks this works on tokens, not on Julia semantics: there
 * is no notion of call arity or scope. A few heuristics keep the obvious false
 * positives out -- field access, assignments, keyword arguments and comments.
 */

export interface JuliaSymbol {
  href: string
  kind: string
  signatures: string[]
  brief: string
}

interface DocstringEntry {
  page: string
  id: string
  kind: string
  signature: string
  brief: string
}

// A docstring as emitted by DocumenterVitepress:
//   <summary><a id='…' href='…'><span class="jlbinding">Module.name</span></a>
//   <Badge … text="Function" /></summary> … </details>
const DOCSTRING_RE =
  /<summary><a id='([^']+)' href='[^']*'><span class="jlbinding">([^<]+)<\/span><\/a>\s*<Badge[^>]*text="([^"]+)"[^>]*\/><\/summary>([\s\S]*?)<\/details>/g

const SIGNATURE_RE = /```julia\r?\n([\s\S]*?)\r?\n```/

const IDENTIFIER_RE = /[A-Za-z_][A-Za-z0-9_]*!?/g

// Qualifiers under which a documented name still refers to the package.
const PACKAGE_QUALIFIERS = new Set(['BeamletOptics', 'BMO'])

// Never link the package name itself.
const EXCLUDED_NAMES = new Set(['BeamletOptics'])

// Marker attribute handed from the `tokens` hook to the `span` hook.
const REF_ATTR = 'data-jl-ref'

const SKIPPED_DIRS = new Set(['node_modules', 'public', 'components', 'assets'])

function listMarkdown(dir: string): string[] {
  const files: string[] = []
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    if (entry.name.startsWith('.') || SKIPPED_DIRS.has(entry.name)) continue
    const full = path.join(dir, entry.name)
    if (entry.isDirectory()) files.push(...listMarkdown(full))
    else if (entry.name.endsWith('.md')) files.push(full)
  }
  return files
}

function plainText(markdown: string): string {
  return markdown
    .replace(/<[^>]+>/g, '')
    .replace(/\[([^\]]*)\]\([^)]*\)/g, '$1')
    .replace(/`([^`]*)`/g, '$1')
    .replace(/\*\*|__/g, '')
    .replace(/\s+/g, ' ')
    .trim()
}

function firstSentence(text: string, max = 200): string {
  const match = text.match(/^(.+?[.!?])(?=\s+[A-Z(`]|$)/)
  const sentence = match ? match[1] : text
  return sentence.length > max ? sentence.slice(0, max - 1).trimEnd() + '…' : sentence
}

// The first prose paragraph after the signature block, cut to its first sentence.
function briefFrom(body: string): string {
  const signature = body.match(SIGNATURE_RE)
  const rest = signature ? body.slice((signature.index ?? 0) + signature[0].length) : body
  for (const chunk of rest.split(/\r?\n\s*\r?\n/)) {
    // A list that directly follows the prose without a blank line belongs to the
    // same chunk; keep only the lines before it.
    const text = chunk.split(/\r?\n\s*(?:[-*+]|\d+\.)\s/)[0].trim().replace(/:$/, ' …')
    if (!text || /^(```|<|\*\*|-|#|:::)/.test(text)) continue
    return firstSentence(plainText(text))
  }
  return ''
}

/**
 * Build the symbol index from the emitted markdown under `mdRoot`.
 *
 * `base` is the VitePress base path, `htmlExt` is `''` for clean URLs and
 * `'.html'` otherwise. Both are applied here because the links are generated as
 * raw HAST and therefore bypass VitePress' own link normalisation.
 */
export function buildJuliaSymbolIndex(
  mdRoot: string,
  base: string,
  htmlExt: string,
): Map<string, JuliaSymbol> {
  const found = new Map<string, DocstringEntry[]>()
  for (const file of listMarkdown(mdRoot)) {
    const page = path.relative(mdRoot, file).replace(/\\/g, '/').replace(/\.md$/, '')
    const text = fs.readFileSync(file, 'utf8')
    for (const [, id, binding, kind, body] of text.matchAll(DOCSTRING_RE)) {
      const name = binding.slice(binding.lastIndexOf('.') + 1)
      if (name.length < 3 || EXCLUDED_NAMES.has(name)) continue
      const entries = found.get(name) ?? []
      entries.push({
        page,
        id,
        kind,
        signature: body.match(SIGNATURE_RE)?.[1].trim() ?? '',
        brief: briefFrom(body),
      })
      found.set(name, entries)
    }
  }

  const index = new Map<string, JuliaSymbol>()
  for (const [name, entries] of found) {
    // Prefer the topical page over the catch-all reference page, as `@ref` does.
    const target = entries.find((e) => e.page !== 'reference') ?? entries[0]
    // The id is taken verbatim from the target page, so it matches the element id
    // even where DocumenterVitepress leaves a Windows path separator in it.
    const fragment = target.id.replace(/ /g, '%20').replace(/"/g, '%22')
    index.set(name, {
      href: `${base}${target.page}${htmlExt}#${fragment}`,
      kind: target.kind,
      signatures: [...new Set(entries.map((e) => e.signature).filter(Boolean))].slice(0, 4),
      brief: entries.map((e) => e.brief).find(Boolean) ?? '',
    })
  }
  return index
}

interface LinkRange {
  start: number
  end: number
  name: string
}

function linkRanges(line: string, index: Map<string, JuliaSymbol>): LinkRange[] {
  const comment = line.indexOf('#')
  const ranges: LinkRange[] = []
  for (const match of line.matchAll(IDENTIFIER_RE)) {
    const name = match[0]
    const start = match.index ?? 0
    const end = start + name.length
    if (!index.has(name)) continue
    if (comment !== -1 && start > comment) continue
    // Field access such as `lens.radius`, unless qualified with the package.
    if (line[start - 1] === '.') {
      const qualifier = line.slice(0, start - 1).match(/([A-Za-z_][A-Za-z0-9_]*)$/)
      if (!qualifier || !PACKAGE_QUALIFIERS.has(qualifier[1])) continue
    }
    // Assignments and keyword arguments (`radius = 5`), but not comparisons (`==`).
    const following = line.slice(end).match(/^\s*(=+)/)
    if (following && following[1] === '=') continue
    ranges.push({ start, end, name })
  }
  return ranges
}

// Split tokens at identifier boundaries so every linkable name is its own token.
// Shiki often lumps identifiers into one token, e.g. `(system, beam)`.
function splitLine(tokens: ThemedToken[], index: Map<string, JuliaSymbol>): ThemedToken[] {
  const text = tokens.map((t) => t.content).join('')
  const ranges = linkRanges(text, index)
  if (ranges.length === 0) return tokens

  const result: ThemedToken[] = []
  let column = 0
  for (const token of tokens) {
    const start = column
    const end = column + token.content.length
    column = end
    // Names spanning several tokens are rare and simply left unlinked.
    const inside = ranges.filter((r) => r.start >= start && r.end <= end)
    if (inside.length === 0) {
      result.push(token)
      continue
    }
    const piece = (from: number, to: number, name?: string): ThemedToken => ({
      ...token,
      content: text.slice(from, to),
      offset: token.offset + (from - start),
      ...(name ? { htmlAttrs: { ...token.htmlAttrs, [REF_ATTR]: name } } : {}),
    })
    let cursor = start
    for (const range of inside) {
      if (range.start > cursor) result.push(piece(cursor, range.start))
      result.push(piece(range.start, range.end, range.name))
      cursor = range.end
    }
    if (cursor < end) result.push(piece(cursor, end))
  }
  return result
}

export function juliaRefTransformer(index: Map<string, JuliaSymbol>): ShikiTransformer {
  const transformer: ShikiTransformer = {
    name: 'julia-docstring-refs',

    tokens(lines) {
      if (this.options.lang !== 'julia' || index.size === 0) return
      return lines.map((tokens) => splitLine(tokens, index))
    },

    span(node, _line, _col, _lineElement, token) {
      const name = token.htmlAttrs?.[REF_ATTR]
      const symbol = name ? index.get(name) : undefined
      if (!symbol) return
      delete node.properties[REF_ATTR]
      delete node.properties.dataJlRef
      const link: Element = {
        type: 'element',
        tagName: 'a',
        properties: {
          href: symbol.href,
          className: ['jl-ref'],
          dataKind: symbol.kind,
          dataSig: symbol.signatures.join('\n'),
          dataBrief: symbol.brief,
        },
        children: [node],
      }
      return link
    },
  }
  return transformer
}

/**
 * Client script for the hover tooltip, injected through the `head` config.
 *
 * The tooltip is attached to `document.body` rather than to the link: code blocks
 * scroll horizontally (`overflow: auto`), which would clip anything positioned
 * inside them. Event delegation on `document` keeps it working across VitePress'
 * client-side navigation. Text goes in through `textContent` only.
 */
export const juliaRefTooltipScript = `(() => {
  let tip = null;
  const hide = () => { if (tip) tip.remove(); };
  const show = (link) => {
    if (!tip) { tip = document.createElement('div'); tip.className = 'jl-tip'; }
    tip.replaceChildren();
    if (link.dataset.sig) {
      const pre = document.createElement('pre');
      pre.textContent = link.dataset.sig;
      tip.append(pre);
    }
    if (link.dataset.brief) {
      const p = document.createElement('p');
      p.textContent = link.dataset.brief;
      tip.append(p);
    }
    if (!tip.childElementCount) return;
    document.body.append(tip);
    const r = link.getBoundingClientRect();
    const left = Math.max(8, Math.min(r.left, window.innerWidth - tip.offsetWidth - 8));
    const below = r.bottom + 6;
    const top = below + tip.offsetHeight > window.innerHeight ? r.top - tip.offsetHeight - 6 : below;
    tip.style.left = left + 'px';
    tip.style.top = top + 'px';
  };
  document.addEventListener('mouseover', (e) => {
    const link = e.target.closest && e.target.closest('a.jl-ref');
    if (link) show(link);
  });
  document.addEventListener('mouseout', (e) => {
    const link = e.target.closest && e.target.closest('a.jl-ref');
    if (link && !link.contains(e.relatedTarget)) hide();
  });
  document.addEventListener('scroll', hide, true);
})();`
