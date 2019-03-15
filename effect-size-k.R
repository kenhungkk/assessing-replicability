cell2k <- function(...) {
  n <- unlist(list(...))
  1 / sqrt(sum(1 / n))
}

# Columns: study.num, k.O, k.R
df.k <-
  rbind(
    # t-tests and ANOVAS ======================================================
    c(1, sqrt(14), sqrt(29)),
    c(2, sqrt(24), sqrt(24)),
    c(3, sqrt(25), sqrt(32)),
    c(4, cell2k(120, 120), cell2k(136, 134)),
    c(5, sqrt(32), sqrt(48)),
    c(6, sqrt(24), sqrt(32)),
    c(7, sqrt(100), sqrt(15)),
    c(8, cell2k(20, 19), cell2k(13, 20)),
    c(10, cell2k(16, 16), cell2k(16, 15)),
    c(11, sqrt(22), sqrt(30)),
    c(15, sqrt(95), sqrt(242)),
    c(19, cell2k(16, 16), cell2k(11, 10)),
    c(20, cell2k(48, 48), cell2k(47, 61)),
    c(24, sqrt(153), sqrt(49)),
    c(27, sqrt(32), sqrt(71)),
    c(29, sqrt(8), sqrt(15)),
    c(32, sqrt(37), sqrt(38)),
    c(33, sqrt(40), sqrt(40)),
    c(36, cell2k(6, 6, 6, 6), cell2k(6, 6, 6, 6)),
    c(49, cell2k(18, 18), cell2k(39, 49)),
    c(52, cell2k(66.5, 66.5), cell2k(59, 54)),
    c(53, cell2k(16, 15), cell2k(29, 46)),
    c(56, cell2k(25.75, 25.75), cell2k(20, 20)),
    c(58, cell2k(123, 63), cell2k(192, 88)),
    c(61, cell2k(59, 51), cell2k(113, 110)),
    c(63, cell2k(71 / 3, 142 / 3), cell2k(48, 100)),
    c(65, cell2k(25, 20, 58, 52), cell2k(31, 33, 24, 47)),
    c(68, cell2k(40, 40, 40, 40), cell2k(68, 59, 46, 53)),
    c(71, cell2k(116, 259), cell2k(77, 100)),
    c(72, cell2k(65.25, 65.25, 65.25, 65.25), cell2k(64, 61, 62, 64)),
    c(81, cell2k(23.5, 23.5, 23.5, 23.5), cell2k(33, 34, 37, 37)),
    c(87, cell2k(13.75, 13.75, 13.75, 13.75), cell2k(12, 14, 13, 12)),
    c(94, cell2k(13, 15), cell2k(29, 32)),
    c(97, cell2k(18.5, 18.5, 20, 20), cell2k(691, 679, 57, 63)),
    c(106, cell2k(24, 12), cell2k(31, 16)),
    c(107, cell2k(41, 45), cell2k(79, 79)),
    c(110, cell2k(178, 102), cell2k(65, 79)),
    c(111, cell2k(14.75, 14.75, 14.75, 14.75), cell2k(30, 30, 30, 30)),
    c(113, sqrt(125), sqrt(176)),
    c(114, cell2k(16, 16), cell2k(16, 16)),
    c(115, sqrt(32), sqrt(8)),
    c(116, sqrt(173), sqrt(140)),
    c(
      118, cell2k(14, 20, 13, 14, 20, 13, 13, 13),
      cell2k(23, 21, 23, 18, 13, 22, 21, 25)
    ),
    c(121, sqrt(12), sqrt(24)),
    c(122, sqrt(8), sqrt(17)),
    c(124, cell2k(18, 18), cell2k(34, 36)),
    c(127, sqrt(29), sqrt(26)),
    c(133, sqrt(24), sqrt(38)),
    c(136, cell2k(15, 15), cell2k(29, 29)),
    c(145, cell2k(20, 20, 20, 20), cell2k(10, 10, 10, 10)),
    c(146, sqrt(15), sqrt(12)),
    c(148, cell2k(50.5, 50.5, 48.5, 48.5), cell2k(93, 50, 75, 45)),
    c(149, cell2k(50.5, 50.5, 48.5, 48.5), cell2k(107, 70, 83, 58)),
    c(150, sqrt(14), sqrt(19)),
    c(151, cell2k(21, 22), cell2k(68, 58)),
    c(153, sqrt(8), sqrt(8)),
    c(158, cell2k(20, 20), cell2k(46, 49)),
    c(161, cell2k(12, 12, 12, 12), cell2k(12, 13, 11, 12)),
    c(167, sqrt(18), sqrt(22)),
    # Correlations and regressions ============================================
    c(39, 6.926, 7.824),
    c(44, sqrt(71 - 5), sqrt(180 - 5)),
    c(48, sqrt(100 - 9), sqrt(200 - 9)),
    c(93, sqrt(91 - 9), sqrt(76 - 9)),
    c(112, sqrt(11 - 3), sqrt(11 - 3)),
    c(120, sqrt(31 - 3), sqrt(43 - 4)),
    c(134, sqrt(119 - 4), sqrt(238 - 4)),
    c(154, sqrt(70 - 3), sqrt(16 - 3)),
    c(155, sqrt(53 - 3), sqrt(72 - 3))
  )

df.k <- as.data.frame(df.k)
colnames(df.k) <- c("study.num", "k.O", "k.R")