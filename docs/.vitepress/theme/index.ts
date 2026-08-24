// https://vitepress.dev/guide/custom-theme
import { h } from 'vue'
import type { Theme } from 'vitepress'
import DefaultTheme from 'vitepress/theme'
// @ts-expect-error Vite resolves Vue SFCs at runtime; the project lacks the Vue module declaration.
import SidebarToggle from './SidebarToggle.vue'
import './style.css'
import "@catppuccin/vitepress/theme/macchiato/blue.css";

export default {
  extends: DefaultTheme,
  Layout: () => {
    return h(DefaultTheme.Layout, null, {
      // https://vitepress.dev/guide/extending-default-theme#layout-slots
    'layout-top': () => h(SidebarToggle)
    })
  },
  enhanceApp({ app, router, siteData }) {
    // ...
  }
} satisfies Theme
