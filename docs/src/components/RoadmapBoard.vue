<script setup>
import { computed } from 'vue'
// Written by `DocUtils.fetch_milestones` before every docs build (not tracked by git): the
// response of the GitHub API for the milestones of the repository, or `[]` if it failed.
import milestones from './milestones.json'

const MILESTONES_URL = 'https://github.com/JuliaPhysics/BeamletOptics.jl/milestones'

// Status glyphs double as the accessible label so progress never depends on
// colour alone (some readers are colour-blind, and dark/light themes shift
// hue anyway).
const LABELS = { done: 'Done', wip: 'In progress', todo: 'Planned' }
const GLYPHS = { done: '✓', wip: '◐', todo: '○' }

// The status is derived from the milestone, never authored: a closed milestone
// is done, an open one is in progress as soon as one of its issues is closed.
function statusOf(milestone) {
  if (milestone.state === 'closed') {
    return 'done'
  }
  return milestone.closed_issues > 0 ? 'wip' : 'todo'
}

function percentOf(closed, total) {
  return total ? Math.round((100 * closed) / total) : 0
}

// Milestones are versions, so sort their titles numerically ("v0.14.10" after
// "v0.14.9"). Everything is computed from the static JSON import (no refs, no
// lifecycle hooks) so the board is fully present in the server-rendered HTML,
// not injected after hydration.
const cards = computed(() =>
  [...milestones]
    .sort((a, b) => a.title.localeCompare(b.title, undefined, { numeric: true }))
    .map((milestone) => {
      const closed = milestone.closed_issues
      const total = milestone.open_issues + closed
      return {
        id: milestone.number,
        title: milestone.title,
        description: milestone.description,
        url: milestone.html_url,
        status: statusOf(milestone),
        closed,
        total,
        percent: milestone.state === 'closed' ? 100 : percentOf(closed, total),
      }
    })
)

const overall = computed(() => {
  const count = (status) => cards.value.filter((card) => card.status === status).length
  const closed = cards.value.reduce((sum, card) => sum + card.closed, 0)
  const total = cards.value.reduce((sum, card) => sum + card.total, 0)
  return {
    done: count('done'),
    wip: count('wip'),
    todo: count('todo'),
    closed,
    total,
    percent: percentOf(closed, total),
  }
})
</script>

<template>
  <section class="bmo-roadmap">
    <div v-if="cards.length === 0" class="bmo-roadmap-summary">
      The milestones could not be loaded for this build, see the
      <a :href="MILESTONES_URL">milestones on GitHub</a>.
    </div>
    <template v-else>
      <div class="bmo-roadmap-summary">
        <div class="bmo-roadmap-summary-head">
          <span>{{ overall.closed }} of {{ overall.total }} issues closed</span>
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
          <span class="bmo-roadmap-bar-fill done" :style="{ width: overall.percent + '%' }" />
        </div>
        <div class="bmo-roadmap-legend">
          <span class="bmo-roadmap-key done">{{ GLYPHS.done }} {{ overall.done }} released</span>
          <span class="bmo-roadmap-key wip">{{ GLYPHS.wip }} {{ overall.wip }} in progress</span>
          <span class="bmo-roadmap-key todo">{{ GLYPHS.todo }} {{ overall.todo }} planned</span>
          <a class="bmo-roadmap-source" :href="MILESTONES_URL">Milestones on GitHub</a>
        </div>
      </div>

      <ol class="bmo-roadmap-track">
        <li v-for="card in cards" :key="card.id" class="bmo-roadmap-phase">
          <span
            class="bmo-roadmap-marker"
            :class="card.status"
            :aria-label="LABELS[card.status]"
          >{{ GLYPHS[card.status] }}</span>
          <div class="bmo-roadmap-card">
            <div class="bmo-roadmap-card-head">
              <h3><a :href="card.url">{{ card.title }}</a></h3>
              <span class="bmo-roadmap-count">{{ card.closed }} / {{ card.total }} issues</span>
            </div>
            <div
              class="bmo-roadmap-bar"
              role="progressbar"
              :aria-valuenow="card.percent"
              aria-valuemin="0"
              aria-valuemax="100"
              :aria-label="card.title + ' progress'"
            >
              <span class="bmo-roadmap-bar-fill done" :style="{ width: card.percent + '%' }" />
            </div>
            <p v-if="card.description" class="bmo-roadmap-description">{{ card.description }}</p>
          </div>
        </li>
      </ol>
    </template>
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

.bmo-roadmap-source {
  margin-left: auto;
}

.bmo-roadmap-description {
  margin: 8px 0 0;
  font-size: 14px;
  color: var(--vp-c-text-2);
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
