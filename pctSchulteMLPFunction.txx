namespace pct
{

SchulteMLPFunction
::SchulteMLPFunction()
{
  // We operate a change of origin, u0 is always 0
  m_u0=0.;

  // Construct the constant part of R0 and R1 (equations 11 and 14)
  m_R0(0,0) = 1.;
  m_R0(1,0) = 0.;
  m_R0(1,1) = 1.;
  m_R1 = m_R0;

  // Transpose
  m_R0T = m_R0.GetTranspose();
  m_R1T = m_R1.GetTranspose();
}

void
SchulteMLPFunction
::Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut)
{
  m_uOrigin = posIn[2];
  //m_IntForSigmaSqTheta0  = Functor::SchulteMLP::IntegralForSigmaSqTheta ::GetValue(m_u0);
  //m_IntForSigmaSqTTheta0 = Functor::SchulteMLP::IntegralForSigmaSqTTheta::GetValue(m_u0);
  //m_IntForSigmaSqT0      = Functor::SchulteMLP::IntegralForSigmaSqT     ::GetValue(m_u0);

  m_u2 = posOut[2]-m_uOrigin;
  m_IntForSigmaSqTheta2  = Functor::SchulteMLP::IntegralForSigmaSqTheta ::GetValue(m_u2);
  m_IntForSigmaSqTTheta2 = Functor::SchulteMLP::IntegralForSigmaSqTTheta::GetValue(m_u2);
  m_IntForSigmaSqT2      = Functor::SchulteMLP::IntegralForSigmaSqT     ::GetValue(m_u2);

  // Parameters vectors
  m_x0[0] = posIn[0];
  m_x0[1] = std::atan(dirIn[0]);  //dirIn[2] is implicitely 1.
  m_x2[0] = posOut[0];
  m_x2[1] = std::atan(dirOut[0]); //dirOut[2] is implicitely 1.

  m_y0[0] = posIn[1];
  m_y0[1] = std::atan(dirIn[1]);  //dirIn[2] is implicitely 1.
  m_y2[0] = posOut[1];
  m_y2[1] = std::atan(dirOut[1]); //dirOut[2] is implicitely 1.
}

void
SchulteMLPFunction
::Evaluate( const double u, double &x, double&y )
{
#ifdef MLP_TIMING
  m_EvaluateProbe1.Start();
#endif
  const double u1 = u-m_uOrigin;

  // Finish constructing rotation matrices (equations 11 and 14)
  m_R0(0,1) = u1;
  m_R1(0,1) = m_u2-u1;
  m_R0T(1,0) = m_R0(0,1);
  m_R1T(1,0) = m_R1(0,1);

  // Constants used in both integrals
  const double intForSigmaSqTheta1  = Functor::SchulteMLP::IntegralForSigmaSqTheta ::GetValue(u1);
  const double intForSigmaSqTTheta1 = Functor::SchulteMLP::IntegralForSigmaSqTTheta::GetValue(u1);
  const double intForSigmaSqT1      = Functor::SchulteMLP::IntegralForSigmaSqT     ::GetValue(u1);

  // Construct Sigma1 (equations 6-9)
  m_Sigma1(1,1) = intForSigmaSqTheta1/* - m_IntForSigmaSqTheta0*/;
  m_Sigma1(0,1) = u1 * m_Sigma1(1,1) - intForSigmaSqTTheta1/* + m_IntForSigmaSqTTheta0*/;
  m_Sigma1(1,0) = m_Sigma1(0,1);
  m_Sigma1(0,0) = u1 * ( 2*m_Sigma1(0,1) - u1*m_Sigma1(1,1) ) + intForSigmaSqT1/* - m_IntForSigmaSqT0*/;
  if(u1 == m_u0)
    m_Sigma1 *= 0;
  else  
    m_Sigma1 *= Functor::SchulteMLP::ConstantPartOfIntegrals::GetValue(m_u0,u1);

  // Construct Sigma2 (equations 15-18)
  m_Sigma2(1,1) = m_IntForSigmaSqTheta2 - intForSigmaSqTheta1;
  m_Sigma2(0,1) = m_u2 * m_Sigma2(1,1) - m_IntForSigmaSqTTheta2 + intForSigmaSqTTheta1;
  m_Sigma2(1,0) = m_Sigma2(0,1);
  m_Sigma2(0,0) = m_u2 * ( 2*m_Sigma2(0,1) - m_u2*m_Sigma2(1,1) ) + m_IntForSigmaSqT2 - intForSigmaSqT1;
  if(u1 == m_u2)
    m_Sigma2 *=0;
  else  
    m_Sigma2 *= Functor::SchulteMLP::ConstantPartOfIntegrals::GetValue(u1,m_u2);

#ifdef MLP_TIMING
  m_EvaluateProbe1.Stop();
  m_EvaluateProbe2.Start();
#endif

  // x and y, equation 24
  // common calculations
  itk::Matrix<double, 2, 2> Sigma1_R1T = m_Sigma1 * m_R1.GetTranspose();
  itk::Matrix<double, 2, 2> R1_Sigma1 = m_R1 * m_Sigma1;

  InverseMatrix(m_R1);
  itk::Matrix<double, 2, 2> part1 = m_R1 * m_Sigma2 + Sigma1_R1T;
  InverseMatrix(part1);
  itk::Matrix<double, 2, 2> part2 = R1_Sigma1 + m_Sigma2 * m_R1.GetTranspose();
  InverseMatrix(part2);

  // x
  itk::Vector<double, 2> xMLP;
  xMLP = m_R1 * m_Sigma2 * part1 * m_R0 * m_x0 + m_Sigma1 * part2 * m_x2;
  x = xMLP[0];

  // y
  itk::Vector<double, 2> yMLP;
  yMLP = m_R1 * m_Sigma2 * part1 * m_R0 * m_y0 + m_Sigma1 * part2 * m_y2;
  y = yMLP[0];

#ifdef MLP_TIMING
  m_EvaluateProbe2.Stop();
#endif
}

void
SchulteMLPFunction
::EvaluateError( const double u, itk::Matrix<double, 2, 2> &error )
{
  const double u1 = u - m_uOrigin;
  if(u1 > m_u0 && u1 < m_u2)
    {
    double x, y;
    Evaluate(u,x,y);
    itk::Matrix<double, 2, 2> Sigma2_InvR1T = m_Sigma2 * m_R1.GetTranspose();
    InverseMatrix(m_R1);
    itk::Matrix<double, 2, 2> part = Sigma2_InvR1T + m_R1 * m_Sigma1;
    InverseMatrix(part);
    InverseMatrix(m_R1);
    error = m_Sigma1 * part * m_Sigma2 * m_R1.GetTranspose();
    }
}

#ifdef MLP_TIMING
void
SchulteMLPFunction
::PrintTiming(std::ostream& os)
{
  os << "SchulteMLPFunction timing:" << std::endl;
  os << "  EvaluateProbe1: " << m_EvaluateProbe1.GetTotal()
     << ' ' << m_EvaluateProbe1.GetUnit() << std::endl;
  os << "  EvaluateProbe2: " << m_EvaluateProbe2.GetTotal()
     << ' ' << m_EvaluateProbe2.GetUnit() << std::endl;
}
#endif

void
SchulteMLPFunction
::InverseMatrix(itk::Matrix<double, 2, 2> &mat)
{
  double det = 1. / ( mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0) );
  std::swap( mat(0,0), mat(1,1) );
  mat(1,0) *= -1.;
  mat(0,1) *= -1.;
  mat *= det;
}

void
SchulteMLPFunction
::SetTrackerInfo(const double dentry, const double dexit, const double dT, const double sigmap, const double xOverX0 )
{
  m_dentry = dentry* CLHEP::mm;
  m_dexit = dexit * CLHEP::mm;
  m_dT = dT * CLHEP::cm;
  m_sigmap = sigmap * CLHEP::mm;
  m_xOverX0 = xOverX0;
}

void
SchulteMLPFunction
::EvaluateErrorWithTrackerUncertainty(const double u, const double eIn, const double eOut, itk::Matrix<double, 2, 2> &error)
{
  double u1 = u - m_uOrigin;
  itk::Matrix<double, 2, 2> SD;
  SD(0,0) = 1;
  SD(0,1) = 0;
  SD(1,0) = 0;
  SD(1,1) = 1;

  if(u1 > m_u2)
    {
    SD(0,1) = std::abs(u1 - m_u2);
    u1 = m_u2;
    }
  else if(u1 < m_u0) 
    {
    SD(0,1) = std::abs(m_u0 - u1);
    u1 = m_u0;
    }

  double proton_mass_c2 = 938.272013 * CLHEP::MeV;
  double x, y;
  Evaluate(u1 + m_uOrigin,x,y);

  itk::Matrix<double, 2, 2> Sin;
  itk::Matrix<double, 2, 2> Sout;
  itk::Matrix<double, 2, 2> SigmaIn;
  itk::Matrix<double, 2, 2> SigmaOut;

  double betapIn = (eIn + 2*proton_mass_c2)* eIn / (eIn + proton_mass_c2);
  double betapOut = (eOut + 2*proton_mass_c2)* eOut / (eOut + proton_mass_c2);

  double sigmascIn = 13.6*CLHEP::MeV / betapIn * std::sqrt(m_xOverX0) * (1 + 0.038 * std::log(m_xOverX0)); 
  double sigmascOut = 13.6*CLHEP::MeV / betapOut * std::sqrt(m_xOverX0) * (1 + 0.038 * std::log(m_xOverX0)); 
  Sin(0,0) = 1;
  Sin(1,1) = 1;
  Sin(1,0) = 0;
  Sout = Sin;
  Sin(0,1) = m_dentry;
  Sout(0,1) = m_dexit;
  SigmaIn(1,0) = -1./m_dT;
  SigmaIn(1,1) = 1./m_dT;
  SigmaOut = SigmaIn;
  SigmaIn(0,0) = 0;
  SigmaIn(0,1) = 1;
  SigmaOut(0,0) = 1;
  SigmaOut(0,1) = 0;
  SigmaIn = SigmaIn * SigmaIn.GetTranspose() * m_sigmap * m_sigmap;
  SigmaIn(1,1) += sigmascIn * sigmascIn;  
  SigmaOut = SigmaOut * SigmaOut.GetTranspose() * m_sigmap * m_sigmap;
  SigmaOut(1,1) += sigmascOut * sigmascOut;

  InverseMatrix(Sout);
  itk::Matrix<double, 2, 2> C1 = m_R0 * Sin * SigmaIn * Sin.GetTranspose() * m_R0T + m_Sigma1;
  itk::Matrix<double, 2, 2> C2 = m_R1 * Sout * SigmaOut * Sout.GetTranspose() * m_R1.GetTranspose() + m_R1 * m_Sigma2 * m_R1.GetTranspose();

  itk::Matrix<double, 2, 2> C1C2= C1 + C2;
  InverseMatrix(C1C2);
  error = C1 * C1C2 * C2;
  if(u - m_uOrigin > m_u2)
    {
    error = SD * error * SD.GetTranspose();
    }
  else if(u - m_uOrigin < m_u0) 
    {
    InverseMatrix(SD);
    error = SD * error * SD.GetTranspose();
    }

}

}
