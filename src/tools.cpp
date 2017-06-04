#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	int est_dim = estimations.size();
	int truth_dim = ground_truth.size();

	if (est_dim != truth_dim || est_dim == 0) {
		throw std::invalid_argument( "Vector dimensions don't match, or incorrect estimation dimension");
	}

	VectorXd rmse = VectorXd::Zero(estimations[0].size());
	
	// Accumulate squared residuals
	VectorXd res;
	VectorXd sq_res;

	for (int i = 0; i < est_dim; ++i ){
	    res = estimations[i] - ground_truth[i];
	    sq_res = res.array() * res.array();
	    rmse = rmse + sq_res;
//		cout << res;
	}

	// Calculate the mean
	
	rmse = rmse / est_dim;
	
	// Calculate the squared root
    rmse = rmse.array().sqrt();
	// std::cout << rmse << std::endl;

	return rmse;
}