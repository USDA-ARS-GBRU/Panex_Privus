import { defineConfig } from "vite";
import { viteSingleFile } from "vite-plugin-singlefile";

// Build a single self-contained index.html (all JS/CSS inlined) so Privy can ship
// it as a package asset and inject a JSON data blob at runtime — no server, no
// Node required by end users. See scratch/notes/51_frontend_stack_modelC.md.
export default defineConfig({
  plugins: [viteSingleFile()],
  build: {
    target: "es2019",
    cssCodeSplit: false,
    assetsInlineLimit: 100000000,
  },
});
