# V2 track reconstruction

The two plotting input tables were reconstructed from the supplied original V2 SVG. Integer coordinates are inferred from arc angles; the maximum rounding residual is below 0.03 bp. Colors and label names are taken from the SVG. These values reproduce the original graphical intervals, not corrected gene boundaries.

The pasted highlight example is from NRHc5028 and is used as a format reference only.

| Feature | Recovered plotted interval | MFannot feature interval |
|---|---|---|
| nad1 | 1–2147 | 20–2151 |
| trnI(gau) | 2152–2234 | 2176–2247 |
| cob1 | 2308–5654 | 2335–5661 |
| trnA(ugc) | 5662–5768 | 5710–5781 |
| trnF(gaa) | 5782–5860 | 5802–5873 |
| trnY(gua) | 5874–5938 | 5880–5963 |
| trnH(gug) | 5964–6028 | 5970–6041 |
| rnpB | 6102–6286 | 6106–6298 |
| trnM(cau)_1 | 6299–6383 | 6325–6397 |
| rps3 | 6458–7124 | 6464–7171 |
| trnS(uga) | 7412–7478 | 7420–7504 |
| rnl | 7685–11382 | Not annotated in the supplied masterfile |
| cox3 | 11649–12502 | 11662–12507 |
| nad4 | 12688–14175 | 12735–14186 |
| trnD(guc) | 14187–14300 | 14242–14314 |
| nad4L | 14435–14756 | 14490–14756 |
| nad5 | 14730–16737 | 14756–16744 |
| trnP(ugg) | 16745–16857 | 16799–16871 |
| trnN(guu) | 16932–17003 | 16945–17016 |
| nad6 | 17077–17683 | 17083–17691 |
| trnM(cau)_2 | 17812–17883 | 17825–17896 |
| trnE(uuc) | 17897–17963 | 17905–17976 |
| trnT(ugu) | 17963–18035 | 17977–18047 |
| trnQ(uug) | 18048–18140 | 18082–18154 |
| trnK(cuu) | 18155–18258 | 18200–18271 |
| trnS(gcu) | 18272–18351 | 18293–18377 |
| trnW(cca) | 18438–18522 | 18464–18535 |
| rns | 18956–20258 | 18998–20306 |
| atp6 | 20367–21120 | 20400–21164 |
| atp9 | 21225–21455 | 21275–21493 |
| cox1 | 21614–24296 | 21674–24314 |
| atp8 | 24375–24539 | 24419–24565 |
| trnL(uag) | 24626–24723 | 24665–24746 |
| trnR(ucg) | 24747–24844 | 24786–24857 |
| trnG(ucc) | 24858–24957 | 24899–24970 |
| nad2 | 24971–26520 | 25021–26520 |
| nad3 | 26461–26881 | 26520–26900 |
| trnV(uac) | 26961–27025 | 26967–27038 |
| cox2 | 27099–27874 | 27154–27909 |

Many recovered boundaries match the start of the preceding sequence line in the masterfile. This is consistent with a previous parsing offset. No correction is applied to the reproduction tables; use the original annotation to review biological feature boundaries.

Rendering validation: Circos 0.69-8 on DCC generated both PNG and SVG successfully. All 147 ideogram/tick elements, all 15 link elements, and all 39 highlight paths (plus their axis group) match the reference SVG exactly. All 18 label names match; some label positions differ because their original input intervals were not retained. See validation/svg_comparison.json for per-label displacements. The original SVG may also contain manual label placement; this has not been established. The rendering log contains non-fatal uninitialized-value warnings from Circos; both outputs were generated and structurally checked.
