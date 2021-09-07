#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

// Space Prolongation operator starts here
Eigen::VectorXf Space_Prolongation(Eigen::VectorXf d, int i){

    Eigen::VectorXf d_space_prolongated;
    
    int n = d.size()/((i-1)*2);

    //std::cout << "\n n: \n" << n;

    Eigen::MatrixXf space_prolongation = Eigen::MatrixXf::Zero( (2*(i-1) + 1) * 2 * n, d.size());

    //std::cout << "\n space_prolongation size: \n" << space_prolongation.rows();
        
    Eigen::MatrixXf d_ = Eigen::MatrixXf::Zero(i-1, 2*(i-1) + 1);
    
    Eigen::MatrixXf d_m(1, 3);

    d_m << 1, 2, 1;
    
    for(int j=0; j < i-1; ++j){

        d_.block(j, 2*j, 1, 3) = d_m;
    }

    //std::cout << "\n d_: \n" << d_;

    Eigen::MatrixXf time_block = Eigen::MatrixXf::Zero((i-1) * 2, (2*(i-1) + 1)* 2);

    time_block.block(0, 0, i-1, 2*(i-1) + 1) = d_;
    time_block.block(i-1, 2*(i-1) + 1, i-1, 2*(i-1) + 1) = d_;

    Eigen::MatrixXf Prolongation_operator = time_block.transpose();

    std::cout << "\n Prolongation_operator: \n" << Prolongation_operator;

    for(int j=0; j < n; ++j){

        space_prolongation.block((2*(i-1) + 1)* 2 * j, (i-1) * 2 * j, (2*(i-1) + 1)* 2, (i-1) * 2) = Prolongation_operator;
    }

    // std::cout << "\n space_prolongation: \n" << space_prolongation;

    d_space_prolongated = space_prolongation * d;

    return d_space_prolongated;
}
// Restriction operator ends here



// Time prolongation operator starts here
Eigen::VectorXf Time_Prolongation(Eigen::VectorXf d, int i){

    Eigen::VectorXf d_time_prolongated;

    float n = d.size()/((i-1)*2);

    Eigen::MatrixXf time_prol_operator = Eigen::MatrixXf::Zero((n*2 - 1)*(i-1)*2, d.size());

    Eigen::MatrixXf time_block = Eigen::MatrixXf::Identity((i-1)*2, (i-1)*2);

    Eigen::MatrixXf time_block_even = Eigen::MatrixXf::Zero((i-1)*2, (i-1)*2*2);

    time_block_even.block(0, 0, (i-1)*2, (i-1)*2) = time_block;
    time_block_even.block(0, (i-1)*2, (i-1)*2, (i-1)*2) = time_block;

    for(int j=0; j < (n * 2 - 1); ++j){

        if(j%2 == 0){
            time_prol_operator.block((i-1)*2*j, (i-1)*j, (i-1)*2, (i-1)*2) = time_block;
        }
        
    }

    for(int j=0; j < (n * 2 - 1) - 1; ++j){

        if(j%2 == 0){
            time_prol_operator.block((i-1)*2*(j+1), (i-1)*j, (i-1)*2, (i-1)*2*2) = (1.0/2) * time_block_even;
        }
        
    }

    d_time_prolongated = (1.0/2) * time_prol_operator * d;

    return d_time_prolongated;
}

// Prolongation starts here
Eigen::VectorXf Prolongation(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf w_prolongated;

    Eigen::VectorXf Time_w_pro = Time_Prolongation(w_old, i);
    //std::cout << "\n time_w_prolongated size: \n" << Time_w_pro.size();

    Eigen::VectorXf Space_w_pro = Space_Prolongation(Time_w_pro, i);
    std::cout << "\n w_res size: \n" << Space_w_pro.size();

    return w_prolongated;
}

int main(){

    int n;
    int i;
    
    std::cout << "Enter the number of time step: " << std::endl;

    std::cin >> n;

    std::cout << "Enter the number of space grids: " << std::endl;

    std::cin >> i;

    float delt = 1.0/n;
    float delx = 1.0/i;

    std::cout << "\n delt: " << delt << std::endl;
    std::cout << "\n delx: " << delx << std::endl;

    Eigen::VectorXf w_old = Eigen::VectorXf::Random((i-1)*2*n + 2*(i-1));

    std::cout << "\n w_old size: \n" << w_old.size() << std::endl;
    //std::cout << "\n w_old: \n" << w_old << std::endl;

    Eigen::VectorXf w_prolongated = Prolongation(w_old, n, i);
}