#ifndef SVG_CHECK_H
#define SVG_CHECK_H

/*
 * Well-formedness checker for the markup of the report plots (test code only; no external XML
 * tool). Shared by svg_plot_test and report_plots_test.
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* well-formedness checker (no external XML tool)                       */
/* ------------------------------------------------------------------ */

#define CHK_DEPTH 64
#define CHK_IDS 256

static int is_name_start(int c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_' || c == ':'; }
static int is_name_char(int c) { return is_name_start(c) || (c >= '0' && c <= '9') || c == '-' || c == '.'; }
static int is_space(int c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }

/* '&' must start one of the entities this library emits (or a numeric reference). */
static int entity_ok(const char* s)
{
        static const char* ent[] = {"&amp;", "&lt;", "&gt;", "&quot;", "&#39;", "&apos;"};
        for (size_t i = 0; i < sizeof(ent) / sizeof(ent[0]); i++) {
                if (strncmp(s, ent[i], strlen(ent[i])) == 0) {
                        return 1;
                }
        }
        return 0;
}

/* Allowed characters of a numeric attribute (non-finite numbers cannot sneak in). */
static int numeric_attr_ok(const char* name, const char* v, size_t n)
{
        const char* allowed = NULL;
        if (strcmp(name, "d") == 0) {
                allowed = "MmLlHhVvZz0123456789.,- ";
        } else if (strcmp(name, "x") == 0 || strcmp(name, "y") == 0 || strcmp(name, "width") == 0 ||
                   strcmp(name, "height") == 0 || strcmp(name, "stroke-width") == 0) {
                allowed = "0123456789.-";
        } else if (strcmp(name, "viewBox") == 0) {
                allowed = "0123456789. -";
        } else if (strcmp(name, "transform") == 0) {
                allowed = "0123456789.- ()abeilnorst";
        }
        if (!allowed) {
                return 1;
        }
        int commas = 0;   /* one comma per x,y pair; two in a row would be a decimal comma */
        for (size_t i = 0; i < n; i++) {
                if (!strchr(allowed, v[i])) {
                        return 0;
                }
                if (v[i] == ',') {
                        commas++;
                } else if (v[i] == ' ' || isalpha((unsigned char)v[i])) {
                        commas = 0;
                }
                if (commas > 1) {
                        return 0;
                }
        }
        return n > 0 || strcmp(name, "d") == 0;   /* an empty grid-line path is fine */
}

/*
 * Returns 1 when doc (a sequence of <svg> elements and whitespace/text) is
 * well formed: tags balanced and properly nested, attributes quoted and unique,
 * entities valid, ids unique, numeric attributes finite, no forbidden
 * elements, and every <svg> has role="img" plus a <title> and a <desc> as
 * direct children. On failure err describes the first problem.
 */
static int svg_check_ex(const char* doc, char* err, size_t errn, int allow_html)
{
        char stack[CHK_DEPTH][24];
        int svg_title[CHK_DEPTH], svg_desc[CHK_DEPTH];
        char ids[CHK_IDS][64];
        int depth = 0, nids = 0, svgs = 0, in_svg = 0;
        static const char* forbidden[] = {"script", "circle", "ellipse", "image", "use", "a", "foreignObject", "style"};
        const char* p = doc;

#define BAD(...) do { snprintf(err, errn, __VA_ARGS__); return 0; } while (0)

        while (*p) {
                if (*p != '<') {
                        if (*p == '&' && !entity_ok(p)) {
                                BAD("bad entity at offset %ld", (long)(p - doc));
                        }
                        if (depth == 0 && !is_space(*p)) {
                                BAD("text outside an element at offset %ld", (long)(p - doc));
                        }
                        p++;
                        continue;
                }
                p++;
                if (*p == '!' || *p == '?') {
                        BAD("comment, doctype or processing instruction at offset %ld", (long)(p - doc));
                }
                if (*p == '/') {
                        char name[24];
                        size_t k = 0;
                        p++;
                        while (is_name_char(*p) && k < sizeof(name) - 1) {
                                name[k++] = *p++;
                        }
                        name[k] = 0;
                        if (*p != '>') {
                                BAD("bad closing tag </%s", name);
                        }
                        p++;
                        if (depth == 0 || strcmp(stack[depth - 1], name) != 0) {
                                BAD("closing tag </%s> does not match %s", name, depth ? stack[depth - 1] : "(nothing)");
                        }
                        depth--;
                        if (strcmp(name, "svg") == 0) {
                                in_svg = 0;
                        }
                        continue;
                }
                if (!is_name_start(*p)) {
                        BAD("bad tag start at offset %ld", (long)(p - doc));
                }
                char name[24];
                size_t k = 0;
                while (is_name_char(*p)) {
                        if (k >= sizeof(name) - 1) {
                                BAD("tag name too long");
                        }
                        name[k++] = *p++;
                }
                name[k] = 0;
                for (size_t i = 0; i < sizeof(forbidden) / sizeof(forbidden[0]); i++) {
                        if (strcmp(name, forbidden[i]) == 0) {
                                BAD("forbidden element <%s>", name);
                        }
                }
                int is_svg = strcmp(name, "svg") == 0;
                if (is_svg && in_svg) {
                        BAD("nested <svg>");
                }
                /* the page wraps figures in <p> (and a note may carry <strong>); nothing else */
                int html_el = allow_html && !in_svg && (strcmp(name, "p") == 0 || (depth > 0 && strcmp(name, "strong") == 0));
                if (!is_svg && !in_svg && !html_el) {
                        BAD("<%s> outside an <svg>", name);
                }
                if (depth == 1 && strcmp(name, "title") == 0) {
                        svg_title[0]++;
                }
                if (depth == 1 && strcmp(name, "desc") == 0) {
                        svg_desc[0]++;
                }

                char seen[32][32];
                int nseen = 0, role_img = 0;
                for (;;) {
                        int had_space = 0;
                        while (is_space(*p)) {
                                p++;
                                had_space = 1;
                        }
                        if (*p == '>' || *p == '/') {
                                break;
                        }
                        if (!had_space || !is_name_start(*p)) {
                                BAD("<%s>: bad attribute at offset %ld", name, (long)(p - doc));
                        }
                        char an[32];
                        k = 0;
                        while (is_name_char(*p)) {
                                if (k >= sizeof(an) - 1) {
                                        BAD("attribute name too long");
                                }
                                an[k++] = *p++;
                        }
                        an[k] = 0;
                        for (int i = 0; i < nseen; i++) {
                                if (strcmp(seen[i], an) == 0) {
                                        BAD("<%s>: duplicate attribute %s", name, an);
                                }
                        }
                        if (nseen < 32) {
                                strcpy(seen[nseen++], an);
                        }
                        if (*p++ != '=') {
                                BAD("<%s %s>: attribute without value", name, an);
                        }
                        char q = *p++;
                        if (q != '"' && q != '\'') {
                                BAD("<%s %s>: unquoted attribute value", name, an);
                        }
                        const char* v = p;
                        while (*p && *p != q) {
                                if (*p == '<') {
                                        BAD("<%s %s>: '<' in attribute value", name, an);
                                }
                                if (*p == '&' && !entity_ok(p)) {
                                        BAD("<%s %s>: bad entity", name, an);
                                }
                                p++;
                        }
                        if (*p != q) {
                                BAD("<%s %s>: unterminated attribute value", name, an);
                        }
                        size_t vn = (size_t)(p - v);
                        p++;
                        if (!numeric_attr_ok(an, v, vn)) {
                                BAD("<%s %s>: non-numeric or non-finite value \"%.*s\"", name, an, (int)(vn > 40 ? 40 : vn), v);
                        }
                        if (strcmp(an, "id") == 0) {
                                char idv[64];
                                snprintf(idv, sizeof(idv), "%.*s", (int)(vn > 63 ? 63 : vn), v);
                                for (int i = 0; i < nids; i++) {
                                        if (strcmp(ids[i], idv) == 0) {
                                                BAD("duplicate id %s", idv);
                                        }
                                }
                                if (nids < CHK_IDS) {
                                        strcpy(ids[nids++], idv);
                                }
                        }
                        if (strcmp(an, "role") == 0 && vn == 3 && strncmp(v, "img", 3) == 0) {
                                role_img = 1;
                        }
                }
                if (is_svg) {
                        svgs++;
                        if (!role_img) {
                                BAD("<svg> without role=\"img\"");
                        }
                }
                int self_close = 0;
                if (*p == '/') {
                        self_close = 1;
                        p++;
                }
                if (*p != '>') {
                        BAD("<%s>: expected '>'", name);
                }
                p++;
                if (!self_close) {
                        if (depth >= CHK_DEPTH) {
                                BAD("nesting too deep");
                        }
                        strcpy(stack[depth], name);
                        if (is_svg) {
                                svg_title[0] = 0;
                                svg_desc[0] = 0;
                                in_svg = 1;
                        }
                        depth++;
                } else if (is_svg) {
                        BAD("empty <svg/>");
                }
                /* an <svg> is complete when its closing tag pops it; check title/desc then */
                if (!self_close && is_svg) {
                        const char* end = strstr(p, "</svg>");
                        if (!end) {
                                BAD("<svg> not closed");
                        }
                        /* title and desc must appear before any drawing, as the first two children */
                        if (strncmp(p, "<title", 6) != 0) {
                                BAD("<svg> does not start with <title>");
                        }
                        const char* d = strstr(p, "</title>");
                        if (!d || strncmp(d + 8, "<desc", 5) != 0) {
                                BAD("<svg> has no <desc> after <title>");
                        }
                }
        }
        if (depth != 0) {
                BAD("%d element(s) left open (<%s>)", depth, stack[depth - 1]);
        }
        if (svgs == 0) {
                BAD("no <svg> found");
        }
        return 1;
#undef BAD
}

static inline int svg_check(const char* doc, char* err, size_t errn)
{
        return svg_check_ex(doc, err, errn, 0);
}

#endif
