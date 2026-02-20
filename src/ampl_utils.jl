# Define getters for meta.
for field in fieldnames(NLPModelMeta) ∪ fieldnames(AmplNLPMeta)
  meth = Symbol("get_", field)
  if field ∈ fieldnames(NLPModelMeta)
    @eval begin
      NLPModels.$meth(meta::AmplNLPMeta) = getproperty(meta, $(QuoteNode(field)))
    end
    # do not re-export functions from NLPModels
  else
    @eval $meth(meta::AmplNLPMeta) = getproperty(meta, $(QuoteNode(field)))
    @eval $meth(nlp::AmplModel) = $meth(nlp.meta)
    @eval $meth(nlp::AmplMPECModel) = $meth(nlp.meta)
    @eval export $meth
  end
end
