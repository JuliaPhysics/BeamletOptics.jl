<script setup>
import { ref, computed } from 'vue'
import { withBase } from 'vitepress'
import componentEntries from './catalog.json'
import beamEntries from './beam_catalog.json'

// `catalog` selects the entries and the page folder the tile images are rendered into
// (`catalog_showcase.jl` and `beam_catalog_showcase.jl`).
const props = defineProps({
  catalog: { type: String, default: 'components' },
})

const catalogs = {
  components: { entries: componentEntries, folder: '/basics/components/' },
  beams: { entries: beamEntries, folder: '/basics/beams/' },
}

const entries = catalogs[props.catalog].entries

// Eagerly import every rendered tile PNG so a plain filename in the JSON can be resolved to
// the URL Vite actually serves it at. Entries whose image is missing (e.g. a typo, or a
// placeholder build that skips rendering) fall back to no image instead of breaking the
// build.
const imgs = import.meta.glob(['../basics/components/catalog_*.png', '../basics/beams/catalog_*.png'], {
  eager: true,
  import: 'default',
})

function imageSrc(name) {
  const suffix = catalogs[props.catalog].folder + name
  const match = Object.entries(imgs).find(([path]) => path.endsWith(suffix))
  return match ? match[1] : null
}

const categories = computed(() => {
  const counts = new Map()
  for (const entry of entries) {
    counts.set(entry.category, (counts.get(entry.category) ?? 0) + 1)
  }
  return [
    { name: 'All', count: entries.length },
    ...Array.from(counts, ([name, count]) => ({ name, count })),
  ]
})

// Initial state shows everything so the pre-hydration/SSR output already contains every
// entry -- the category filter is a progressive enhancement, not a requirement to see the
// catalog.
const activeCategory = ref('All')

const filteredEntries = computed(() => {
  if (activeCategory.value === 'All') {
    return entries
  }
  return entries.filter((entry) => entry.category === activeCategory.value)
})
</script>

<template>
  <div class="bmo-catalog-chips">
    <button
      v-for="cat in categories"
      :key="cat.name"
      type="button"
      class="bmo-catalog-chip"
      :class="{ active: activeCategory === cat.name }"
      @click="activeCategory = cat.name"
    >
      {{ cat.name }} <span class="bmo-catalog-chip-count">{{ cat.count }}</span>
    </button>
  </div>
  <div class="bmo-gallery">
    <div v-for="entry in filteredEntries" :key="entry.name" class="bmo-tile">
      <img v-if="imageSrc(entry.image)" :src="imageSrc(entry.image)" :alt="entry.name" />
      <p><a :href="withBase(entry.href)">{{ entry.name }}</a></p>
      <p class="bmo-teaser">{{ entry.teaser }}</p>
      <p class="bmo-tags">{{ entry.category }}</p>
    </div>
  </div>
</template>

<style scoped>
.bmo-catalog-chips {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  margin: 16px 0;
}

.bmo-catalog-chip {
  border: 1px solid var(--vp-c-divider);
  border-radius: 999px;
  background: var(--vp-c-bg-soft);
  color: var(--vp-c-text-1);
  padding: 4px 12px;
  font-size: 13px;
  cursor: pointer;
  transition: border-color 0.2s, color 0.2s;
}

.bmo-catalog-chip:hover {
  border-color: var(--vp-c-brand-1);
}

.bmo-catalog-chip.active {
  border-color: var(--vp-c-brand-1);
  color: var(--vp-c-brand-1);
  font-weight: 600;
}

.bmo-catalog-chip-count {
  color: var(--vp-c-text-3);
}
</style>
