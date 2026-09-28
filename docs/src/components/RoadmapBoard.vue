<script setup>
import { computed, h } from 'vue'
import data from './roadmap.json'

// Status glyphs double as the accessible label so progress never depends on
// colour alone (some readers are colour-blind, and dark/light themes shift
// hue anyway).
const LABELS = { done: 'Done', wip: 'In progress', todo: 'Planned' }
const GLYPHS = { done: '✓', wip: '◐', todo: '○' }

// Progress must be derived from item status rather than hand-maintained, or
// the counters go stale the moment a single item flips -- see roadmap.json.
function leaves(item) {
  return item.children ? item.children.flatMap(leaves) : [item]
}

// A parent's own status is never authored; it always follows its leaves so
// the two can never disagree.
function statusOf(item) {
  if (!item.children) {
    return item.status
  }
  const childLeaves = leaves(item)
  if (childLeaves.every((leaf) => leaf.status === 'done')) {
    return 'done'
  }
  if (childLeaves.every((leaf) => leaf.status === 'todo')) {
    return 'todo'
  }
  return 'wip'
}

// Leaves only: counting parents too would double-count e.g. "Polarizing
// optics" once for itself and once per child.
function tally(items) {
  const allLeaves = items.flatMap(leaves)
  const done = allLeaves.filter((leaf) => leaf.status === 'done').length
  const wip = allLeaves.filter((leaf) => leaf.status === 'wip').length
  const todo = allLeaves.filter((leaf) => leaf.status === 'todo').length
  const total = allLeaves.length
  return { done, wip, todo, total, percent: total ? Math.round((100 * done) / total) : 0 }
}

// Everything below is computed from the static JSON import (no refs, no
// lifecycle hooks) so the board is fully present in the server-rendered
// HTML, not injected after hydration.
const phases = computed(() =>
  data.phases.map((phase) => ({
    ...phase,
    stats: tally(phase.items),
    status: statusOf({ children: phase.items }),
  }))
)

const overall = computed(() => tally(data.phases.flatMap((phase) => phase.items)))

// The nesting in roadmap.json is not fixed at one level (see "Polarizing
// optics"), so the item list is a render-function component that recurses
// into `item.children` instead of a template hard-coded to two levels.
const RoadmapItem = {
  name: 'RoadmapItem',
  props: { item: { type: Object, required: true } },
  setup(props) {
    return () => {
      const { item } = props
      const status = statusOf(item)
      const kids = [
        h(
          'span',
          { class: ['bmo-roadmap-dot', status], 'aria-label': LABELS[status] },
          GLYPHS[status]
        ),
        item.title,
      ]
      if (item.note) {
        kids.push(h('span', { class: 'bmo-roadmap-note' }, item.note))
      }
      if (item.children) {
        kids.push(h(RoadmapItemList, { items: item.children }))
      }
      return h('li', kids)
    }
  },
}

const RoadmapItemList = {
  name: 'RoadmapItemList',
  props: { items: { type: Array, required: true } },
  setup(props) {
    return () =>
      h(
        'ul',
        { class: 'bmo-roadmap-items' },
        props.items.map((item) => h(RoadmapItem, { key: item.title, item }))
      )
  },
}
</script>

<template>
  <section class="bmo-roadmap">
    <div class="bmo-roadmap-summary">
      <div class="bmo-roadmap-summary-head">
        <span>{{ overall.done }} of {{ overall.total }} planned features implemented</span>
        <span class="bmo-roadmap-percent">{{ overall.percent }}%</span>
      </div>
      <div
        class="bmo-roadmap-bar"
        role="progressbar"
        :aria-valuenow="overall.percent"
        aria-valuemin="0"
        aria-valuemax="100"
        aria-label="Overall progress"
      >
        <span
          class="bmo-roadmap-bar-fill done"
          :style="{ width: overall.total ? (100 * overall.done) / overall.total + '%' : '0%' }"
        />
        <span
          class="bmo-roadmap-bar-fill wip"
          :style="{ width: overall.total ? (100 * overall.wip) / overall.total + '%' : '0%' }"
        />
      </div>
      <div class="bmo-roadmap-legend">
        <span class="bmo-roadmap-key done">{{ GLYPHS.done }} {{ overall.done }} done</span>
        <span class="bmo-roadmap-key wip">{{ GLYPHS.wip }} {{ overall.wip }} in progress</span>
        <span class="bmo-roadmap-key todo">{{ GLYPHS.todo }} {{ overall.todo }} planned</span>
      </div>
    </div>

    <ol class="bmo-roadmap-track">
      <li v-for="phase in phases" :key="phase.id" class="bmo-roadmap-phase">
        <span
          class="bmo-roadmap-marker"
          :class="phase.status"
          :aria-label="LABELS[phase.status]"
        >{{ GLYPHS[phase.status] }}</span>
        <div class="bmo-roadmap-card">
          <div class="bmo-roadmap-card-head">
            <h3>{{ phase.title }}</h3>
            <span class="bmo-roadmap-count">{{ phase.stats.done }} / {{ phase.stats.total }}</span>
          </div>
          <div
            class="bmo-roadmap-bar"
            role="progressbar"
            :aria-valuenow="phase.stats.percent"
            aria-valuemin="0"
            aria-valuemax="100"
            :aria-label="phase.title + ' progress'"
          >
            <span
              class="bmo-roadmap-bar-fill done"
              :style="{ width: phase.stats.percent + '%' }"
            />
          </div>
          <RoadmapItemList :items="phase.items" />
        </div>
      </li>
    </ol>
  </section>
</template>

<style scoped>
.bmo-roadmap {
  --bmo-done: var(--julia-green);
  --bmo-wip: var(--julia-purple);
  --bmo-todo: var(--vp-c-text-3);
  --bmo-track: var(--vp-c-divider);
  margin: 16px 0;
}

/* Same lighter shades overrides.css already uses for --vp-c-tip-1 /
   --vp-c-warning-1 in dark mode, so the board matches the custom blocks. */
.dark .bmo-roadmap {
  --bmo-done: #5fb84f;
  --bmo-wip: #b48ecb;
}

.bmo-roadmap-summary {
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  background: var(--vp-c-bg-soft);
  padding: 16px 20px;
  margin-bottom: 24px;
}

.bmo-roadmap-summary-head {
  display: flex;
  justify-content: space-between;
  align-items: baseline;
  font-weight: 600;
  margin-bottom: 8px;
}

.bmo-roadmap-percent {
  color: var(--vp-c-text-2);
}

.bmo-roadmap-bar {
  display: flex;
  height: 8px;
  border-radius: 999px;
  background: var(--bmo-track);
  overflow: hidden;
}

.bmo-roadmap-bar-fill {
  display: block;
  height: 100%;
  transition: width 0.3s;
}

@media (prefers-reduced-motion: reduce) {
  .bmo-roadmap-bar-fill {
    transition: none;
  }
}

.bmo-roadmap-bar-fill.done {
  background: var(--bmo-done);
}

.bmo-roadmap-bar-fill.wip {
  background: var(--bmo-wip);
}

.bmo-roadmap-legend {
  display: flex;
  flex-wrap: wrap;
  gap: 16px;
  margin-top: 10px;
  font-size: 13px;
  color: var(--vp-c-text-2);
}

.bmo-roadmap-key.done {
  color: var(--bmo-done);
}

.bmo-roadmap-key.wip {
  color: var(--bmo-wip);
}

.bmo-roadmap-key.todo {
  color: var(--bmo-todo);
}

/* The rail is a pseudo-element behind the markers rather than a real
   element, so it can span the full track height without affecting layout. */
.bmo-roadmap-track {
  position: relative;
  list-style: none;
  margin: 0;
  padding-left: 40px;
}

.bmo-roadmap-track::before {
  content: '';
  position: absolute;
  left: 11px;
  top: 0;
  bottom: 0;
  width: 2px;
  background: var(--bmo-track);
}

.bmo-roadmap-phase {
  position: relative;
  margin-bottom: 24px;
}

.bmo-roadmap-phase:last-child {
  margin-bottom: 0;
}

.bmo-roadmap-marker {
  position: absolute;
  left: -40px;
  top: 0;
  width: 24px;
  height: 24px;
  display: flex;
  align-items: center;
  justify-content: center;
  border-radius: 50%;
  border: 2px solid var(--bmo-track);
  background: var(--vp-c-bg);
  font-size: 13px;
  line-height: 1;
}

.bmo-roadmap-marker.done {
  color: var(--bmo-done);
  border-color: var(--bmo-done);
}

.bmo-roadmap-marker.wip {
  color: var(--bmo-wip);
  border-color: var(--bmo-wip);
}

.bmo-roadmap-marker.todo {
  color: var(--bmo-todo);
  border-color: var(--bmo-todo);
}

/* Mirrors .bmo-tile in overrides.css so the timeline cards match the rest
   of the docs' card language. */
.bmo-roadmap-card {
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  background: var(--vp-c-bg-soft);
  padding: 16px 20px;
}

.bmo-roadmap-card-head {
  display: flex;
  justify-content: space-between;
  align-items: baseline;
  gap: 8px;
  margin-bottom: 8px;
}

.bmo-roadmap-card-head h3 {
  margin: 0;
  border: none;
  padding: 0;
  font-size: 16px;
}

.bmo-roadmap-count {
  color: var(--vp-c-text-2);
  font-size: 13px;
  white-space: nowrap;
}

/* The item list is built with manual render functions (see script) so it can
   recurse into arbitrarily deep `children`, which means these nodes are never
   compiled from this SFC's own <template> and would not receive the scoped
   `data-v-*` attribute on their own. `:deep()` keeps the styles scoped to
   `.bmo-roadmap` while still reaching that non-template-compiled content. */
.bmo-roadmap :deep(.bmo-roadmap-items) {
  list-style: none;
  margin: 8px 0 0;
  padding-left: 20px;
}

.bmo-roadmap :deep(.bmo-roadmap-items .bmo-roadmap-items) {
  margin-top: 4px;
}

/* Hanging indent: the dot sits in the item's own left padding, so a title that
   wraps (and the note below it) lines up with the text rather than with the
   dot. `text-indent` inherits, so every block that starts a new line inside the
   item has to undo it again. */
.bmo-roadmap :deep(.bmo-roadmap-items li) {
  padding: 2px 0 2px 1.4em;
  text-indent: -1.4em;
}

.bmo-roadmap :deep(.bmo-roadmap-items .bmo-roadmap-items),
.bmo-roadmap :deep(.bmo-roadmap-note) {
  text-indent: 0;
}

.bmo-roadmap :deep(.bmo-roadmap-dot) {
  display: inline-block;
  width: 1.4em;
  font-size: 11px;
  text-align: center;
}

.bmo-roadmap :deep(.bmo-roadmap-dot.done) {
  color: var(--bmo-done);
}

.bmo-roadmap :deep(.bmo-roadmap-dot.wip) {
  color: var(--bmo-wip);
}

.bmo-roadmap :deep(.bmo-roadmap-dot.todo) {
  color: var(--bmo-todo);
}

.bmo-roadmap :deep(.bmo-roadmap-note) {
  display: block;
  font-size: 12px;
  color: var(--vp-c-text-3);
}

@media (max-width: 640px) {
  .bmo-roadmap-track {
    padding-left: 28px;
  }

  .bmo-roadmap-track::before {
    left: 7px;
  }

  .bmo-roadmap-marker {
    left: -28px;
    width: 18px;
    height: 18px;
    font-size: 11px;
  }

  .bmo-roadmap-card-head {
    flex-wrap: wrap;
  }
}
</style>
