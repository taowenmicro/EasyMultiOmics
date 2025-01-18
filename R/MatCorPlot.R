MatCorPlot= function (env.dat, tabOTU, distance = "bray", method = "metal",
          method.cor = "spearman", cor.p = 0.05, x = TRUE, y = TRUE,
          diag = TRUE, sig = TRUE, siglabel = FALSE, shownum = TRUE, numpoint = 22,
          numpoint2 = 21, numsymbol = NULL, curvature = 0.2, lacx = "left",
          lacy = "bottom", range = 0.5, p.thur = 0.3, onlysig = TRUE)
{
  rep = MetalTast(env.dat = env.dat, tabOTU = tabOTU, distance = distance,
                  method = method)
  repR = rep[c(-seq(from = 1, to = dim(rep)[2], by = 2)[-1])]
  repP = rep[seq(from = 1, to = dim(rep)[2], by = 2)]
  ??Miccorplot

  result <- Miccorplot(data = env.dat, method.cor = method.cor,
                       cor.p = cor.p, x = x, y = y, diag = diag, lacx = lacx,
                       lacy = lacy, sig = sig, siglabel = siglabel, shownum = shownum,
                       numpoint = numpoint, numsymbol = numsymbol)

  p = cor_link(data = result[[2]], p = result[[1]], envdata = repR,
           Ptab = repP, numpoint2 = numpoint2, curvature = curvature,
           range = range, lacx = lacx, lacy = lacy, p.thur = p.thur,
           onlysig = onlysig)
  return(p)
}


