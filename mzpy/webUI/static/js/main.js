/* mzpy Web UI – minimal vanilla JS helpers */

(function () {
  'use strict';

  // Auto-hide flash messages after 5 s
  document.addEventListener('DOMContentLoaded', function () {
    const flashes = document.querySelectorAll('.flash');
    flashes.forEach(function (el) {
      setTimeout(function () {
        el.style.opacity = '0';
        el.style.transition = 'opacity 0.5s';
        setTimeout(function () { el.remove(); }, 500);
      }, 5000);
    });
  });

  // Confirm destructive actions
  document.addEventListener('click', function (e) {
    const btn = e.target.closest('[data-confirm]');
    if (!btn) return;
    const msg = btn.getAttribute('data-confirm');
    if (!confirm(msg)) {
      e.preventDefault();
    }
  });
})();
