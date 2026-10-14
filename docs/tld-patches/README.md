# tld patches

`0001-Fix-tld_strbuf_resize-never-updated-alloc_len.patch` fixes
`tld_strbuf_resize` never updating `alloc_len` (every append after the first
growth reallocated the whole buffer). It adds the test `unit-strbuf-growth`
(1e6 appends, `alloc_len >= len`, fewer than 64 reallocations). It was made
against tld commit 961795c and passes tld's full ctest (36 tests).

Maintainer steps (agents never push):

1. `git clone https://github.com/timolassmann/tld && cd tld`
2. `git am /path/to/samstat/docs/tld-patches/0001-*.patch`
3. `cmake -S . -B build && cmake --build build && ctest --test-dir build`
4. Push to tld and note the commit hash.
5. In samstat: `git -C src/tld fetch && git -C src/tld checkout <hash>`, commit
   the submodule pin (this is WP1.5).

WP1.5 then removes the `JSON_BUF_BYTES` chunked flushing from `json_report.c`.
