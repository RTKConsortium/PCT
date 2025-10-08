namespace pct
{

template <class TInputImage, class TOutputImage, class TFFTPrecision>
FDKDDConeBeamVarianceReconstructionFilter<TInputImage, TOutputImage, TFFTPrecision>::
  FDKDDConeBeamVarianceReconstructionFilter()
{
  this->SetNumberOfRequiredInputs(1);

  // Create each filter of the composite filter
  this->m_ExtractFilter = ExtractFilterType::New();
  this->m_WeightFilter = WeightFilterType::New();
  this->m_RampFilter = VarianceFilterType::New();
  this->m_BackProjectionFilter = BackProjectionFilterType::New();

  // Permanent internal connections
  this->m_WeightFilter->SetInput(this->m_ExtractFilter->GetOutput());
  this->m_RampFilter->SetInput(this->m_WeightFilter->GetOutput());
  this->m_BackProjectionFilter->SetProjectionStack(this->m_RampFilter->GetOutput());

  // Default parameters
#if ITK_VERSION_MAJOR >= 4
  this->m_ExtractFilter->SetDirectionCollapseToSubmatrix();
#endif
  this->m_WeightFilter->InPlaceOn();
  this->m_BackProjectionFilter->InPlaceOn();
  this->m_BackProjectionFilter->SetTranspose(true);

  this->m_WeightFilter->SetVarianceReconstruction(true);
}


} // end namespace pct
