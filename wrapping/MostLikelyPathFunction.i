%include "typemaps.i"
%include "std_vector.i"

// pct::MostLikelyPathFunction::Evaluate
%apply double &OUTPUT {double & x, double & y, double & dx, double & dy};
%apply std::vector<double> &OUTPUT {std::vector<double> & x, std::vector<double> & y};