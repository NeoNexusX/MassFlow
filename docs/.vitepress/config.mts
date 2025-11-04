import { defineConfig } from 'vitepress'

/**
 * Site Config with i18n locales
 *
 * This configuration enables bilingual documentation with English as root and Chinese under /zh.
 * To satisfy the current VitePress version typings, per-locale theme configuration
 * is not used; instead, a single themeConfig is defined.
 *
 * Exceptions: None.
 */
export default defineConfig({
  title: 'MassFlow',
  description: 'MassFlow',
  base: '/MassFlow/',
  /**
   * URL and path behavior
   *
   * - cleanUrls: removes trailing `.html` in dev/production for cleaner links.
   * - rewrites: maps English content from `/en/...` to root (`/...`) so English is truly the root locale.
   *
   * Parameters: None.
   * Returns: Static configuration object.
   * Exceptions: None.
   */
  cleanUrls: true,
  rewrites: {
    /** Map English content folder to root paths */
    'en/:path': ':path'
  },
  locales: {
    /** English as the root locale */
    root: {
      label: 'English',
      lang: 'en-US'
    },
    /** Simplified Chinese locale under /zh/ */
    zh: {
      label: '简体中文',
      lang: 'zh-CN'
    }
  },
  markdown: {
    theme: {
      light: 'catppuccin-latte',
      dark: 'catppuccin-macchiato'
    }
  },
  themeConfig: {
    // Single theme config (compatible with current typings)
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Getting Started', link: '/getting-started' },
      { text: 'Contribution', link: '/contribution' }
    ],
    /**
     * Sidebar configuration
     *
     * This lists all currently available documentation pages in both English and Chinese locales.
     *
     * Parameters: None.
     * Returns: Not applicable (static configuration object).
     * Exceptions: None.
     */
    sidebar: [
      {
        text: 'English',
        items: [
          { text: 'Home', link: '/' },
          { text: 'Getting Started', link: '/getting-started' },
          { text: 'Contribution', link: '/contribution' },
          { text: 'Collaboration Guide', link: '/collaboration_guide' },
          { text: 'Naming Conventions', link: '/naming-conventions' }
        ]
      },
    ],
    socialLinks: [
      { icon: 'github', link: 'https://github.com/NeoNexusX/MassFlow' }
    ],
    footer: {
      message: 'Released under the GNU License.',
      copyright: 'Copyright © 2025-present <a href="https://bionet.xmu.edu.cn/">Bionet Team</a>'
    },
    search: {
      provider: 'local'
    }
  }
})
