// View-scoped cleanup for the hash-routed single-page app. Browsers keep ResizeObservers,
// media-query listeners, and timers alive after their DOM node is removed, so a route change
// must dispose them explicitly rather than waiting for a future resize or scroll event.

const cleanupByNode = new WeakMap();
const disposed = new WeakSet();

function run(cleanup) {
  // A teardown failure must never prevent the next route from rendering. Cleanup functions are
  // deliberately small and idempotent; ignoring one failed best-effort teardown is safer than
  // leaving a half-replaced view on screen.
  try { cleanup(); } catch { /* best effort */ }
}

function belongsToDisposedTree(node) {
  for (let current = node; current; current = current.parentNode) {
    if (disposed.has(current)) return true;
  }
  return false;
}

/** Register work that belongs to a DOM node's lifetime. Returns an unregister function. */
export function onDispose(node, cleanup) {
  if (!node || typeof cleanup !== 'function') throw new TypeError('onDispose needs a node and function');
  if (belongsToDisposedTree(node)) {
    run(cleanup);
    return () => {};
  }
  let cleanups = cleanupByNode.get(node);
  if (!cleanups) {
    cleanups = new Set();
    cleanupByNode.set(node, cleanups);
  }
  cleanups.add(cleanup);
  return () => cleanups.delete(cleanup);
}

/** Dispose registered resources below a root, deepest node first. */
export function disposeTree(root, { includeRoot = true } = {}) {
  if (!root) return;
  const nodes = includeRoot ? [root] : [];
  if (root.querySelectorAll) nodes.push(...root.querySelectorAll('*'));
  for (let i = nodes.length - 1; i >= 0; i--) {
    const node = nodes[i];
    disposed.add(node);
    const cleanups = cleanupByNode.get(node);
    if (!cleanups) continue;
    cleanupByNode.delete(node);
    [...cleanups].reverse().forEach(run);
  }
}
