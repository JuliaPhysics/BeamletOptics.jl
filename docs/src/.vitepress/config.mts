import { defineConfig } from 'vitepress'
import { withMermaid } from 'vitepress-plugin-mermaid'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import { buildJuliaSymbolIndex, juliaRefTransformer, juliaRefTooltipScript } from './julia-ref-transformer'
import footnote from "markdown-it-footnote";
import path from 'path'

const mathjax = mathjaxPlugin()

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
}

// Clean URLs only on the deployed site: local static servers such as the VS Code
// Live Server cannot resolve extensionless paths, so a reload there would 404.
const cleanUrls = process.env.CI === 'true'

// Docstring index for the reference links and tooltips in Julia code blocks. It is
// read from the markdown DocumenterVitepress has emitted next to this folder.
const juliaSymbols = buildJuliaSymbolIndex(
  path.resolve(__dirname, '..'),
  baseTemp.base,
  cleanUrls ? '' : '.html',
)

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

// DocumenterVitepress builds the navbar from `pages` with the same nesting as the
// sidebar. A navbar dropdown only holds links and one level of groups of links, so
// deeper groups (e.g. "Geometry" in "API design") are inlined into their parent
// group here. The sidebar keeps the full hierarchy.
function flattenNavGroups(items: any[], depth = 0): any[] {
  return items.flatMap((item) => {
    if (!item.items) return [item]
    const children = flattenNavGroups(item.items, depth + 1)
    return depth >= 2 ? children : [{ ...item, items: children }]
  })
}

const nav = [
  ...flattenNavGroups(navTemp.nav as any),
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default withMermaid(defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}],
    ['script', {}, juliaRefTooltipScript],
    // REPLACE_ME_DOCUMENTER_VITEPRESS_NOINDEX
  ],

  markdown: {
    lineNumbers: true,
    codeTransformers: [juliaReplTransformer(), juliaRefTransformer(juliaSymbols)],
    config(md) {
      md.use(tabsMarkdownPlugin);
      md.use(footnote);
      mathjax.markdownConfig(md);
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  vite: {
    plugins: [
      mathjax.vitePlugin,
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  themeConfig: {
    outline: 'deep',
    logo: { light: '/logo.svg', dark: '/logo-dark.svg' },
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    sidebarDrawer: 'REPLACE_ME_DOCUMENTER_VITEPRESS_SIDEBAR_DRAWER',
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
}))
