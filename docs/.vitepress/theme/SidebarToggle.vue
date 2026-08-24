<script setup lang="ts">
import { onMounted, ref } from 'vue'

const collapsed = ref(false)

function applySidebarState() {
  document.documentElement.classList.toggle(
    'sidebar-collapsed',
    collapsed.value
  )
}

function toggleSidebar() {
  collapsed.value = !collapsed.value

  applySidebarState()

  localStorage.setItem(
    'massflow-sidebar-collapsed',
    String(collapsed.value)
  )
}

onMounted(() => {
  collapsed.value =
    localStorage.getItem('massflow-sidebar-collapsed') === 'true'

  applySidebarState()
})
</script>

<template>
  <button
    class="sidebar-toggle"
    type="button"
    :aria-label="
      collapsed ? 'Expand sidebar' : 'Collapse sidebar'
    "
    :title="
      collapsed ? 'Expand sidebar' : 'Collapse sidebar'
    "
    @click="toggleSidebar"
  >
    {{ collapsed ? '›' : '‹' }}
  </button>
</template>