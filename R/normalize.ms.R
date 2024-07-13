

# ps.ms2 = normalize.ms(ps = ps.ms,method = "rela")
normalize.ms = function(ps,
                        method = "rela"
){
  if (method == "rela") {
    ps1  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x,na.rm =TRUE) )

  }

  if (method == "sampling" ) {
    total = mean(sample_sums(ps));total
    standf = function(x,t = total)round(t*(x/sum(x)))
    ps1 = phyloseq::transform_sample_counts(ps,standf)
  }
  if (method == "log") {
    ps1  = phyloseq::transform_sample_counts(ps, function(x) log(1 + x) )

  }
  return(ps1)
}
