#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <tuple>

// Space Restriction operator starts here
Eigen::MatrixXf Restriction_Operator(int n, int i){

    Eigen::VectorXf d_space_restricted;
        
    Eigen::MatrixXf d_ = Eigen::MatrixXf::Zero((i/2 - 1), i-1);
    
    Eigen::MatrixXf d_m(1, 3);

    d_m << 1, 2, 1;
    
    for(int j=0; j < (i/2 - 1); ++j){

        d_.block(j, 2*j, 1, 3) = d_m;
    }

    Eigen::MatrixXf Operator = Eigen::MatrixXf::Zero((i/2 - 1) * 2, (i-1) * 2);

    Operator.block(0, 0, (i/2 - 1), i-1) = d_;
    Operator.block((i/2 - 1), i-1, (i/2 - 1), i-1) = d_;

    return Operator;
}
// Space Restriction operator ends here

// Space Restriction operator starts here
Eigen::VectorXf Space_Restriction(Eigen::VectorXf d, int i){

    Eigen::VectorXf d_space_restricted;
    
    int n = d.size()/((i-1)*2);

    Eigen::MatrixXf space_res_operator = Eigen::MatrixXf::Zero( (i/2 - 1) * 2 * n, d.size());
    
    Eigen::MatrixXf d_ = Eigen::MatrixXf::Zero((i/2 - 1), i-1);
    
    Eigen::MatrixXf d_m(1, 3);

    d_m << 1, 2, 1;
    
    for(int j=0; j < (i/2 - 1); ++j){

        d_.block(j, 2*j, 1, 3) = d_m;
    }

    Eigen::MatrixXf time_block = Eigen::MatrixXf::Zero((i/2 - 1) * 2, (i-1) * 2);

    time_block.block(0, 0, (i/2 - 1), i-1) = d_;
    time_block.block((i/2 - 1), i-1, (i/2 - 1), i-1) = d_;

    std::cout << "\n restriction operator: \n" << time_block;

    for(int j=0; j < n; ++j){

        space_res_operator.block((i/2 - 1) * 2 * j, (i-1) * 2 * j, (i/2 - 1) * 2, (i-1) * 2) = time_block;
    }

    //std::cout << "\n space_res_operator: \n" << space_res_operator;

    d_space_restricted = space_res_operator * d;

    return d_space_restricted;
}
// Restriction operator ends here

// Time Restriction operator starts here
Eigen::VectorXf Time_Restriction(Eigen::VectorXf d, int n, int i){

    Eigen::VectorXf d_time_restricted;
      
    Eigen::MatrixXf time_res_operator = Eigen::MatrixXf::Zero(((n+1)/2 + 1) * 2 * (i-1), d.size());
    
    Eigen::MatrixXf a_0 = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2 * 2);

    a_0.block(0, 0, (i-1) * 2, (i-1) * 2) = 2 * Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);
    a_0.block(0, (i-1) * 2, (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);

    Eigen::MatrixXf a_m = Eigen::MatrixXf::Zero((i-1) * 2, (i-1) * 2 * 3);

    a_m.block(0, 0, (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);
    a_m.block(0, (i-1) * 2, (i-1) * 2, (i-1) * 2) = 2 * Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);
    a_m.block(0, (i-1) * 2 * 2, (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);

    Eigen::MatrixXf a_nt = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2 * 2);

    a_nt.block(0, 0, (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);
    a_nt.block(0, (i-1) * 2, (i-1) * 2, (i-1) * 2) = 2 * Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);

    time_res_operator.block(0, 0, (i-1) * 2, (i-1) * 2 * 2) = a_0;

    if(n%2 == 0){
        time_res_operator.block(time_res_operator.rows() - (i-1) * 2, time_res_operator.cols() - ((i-1) * 2 * 2), (i-1) * 2, (i-1) * 2 * 2) = a_nt;

        for(int j=1; j < n/2; ++j){

        time_res_operator.block(2*j*(i-1), 2*j*(i-1) + 6*(j-1), (i-1) * 2, (i-1) * 2 * 3) = a_m;

    }    
    }
    else{
        time_res_operator.block(time_res_operator.rows() - (i-1) * 2, time_res_operator.cols() - ((i-1) * 2), (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);;
        
        for(int j=1; j < n/2 + 1; ++j){

        time_res_operator.block(2*j*(i-1), 2*j*(i-1) + 6*(j-1), (i-1) * 2, (i-1) * 2 * 3) = a_m;

    }
    }
    
    time_res_operator = (1.0/4) * time_res_operator;

    //std::cout << "\n time res operator: \n" << time_res_operator;

    d_time_restricted = time_res_operator * d;

    return d_time_restricted;
}
// Restriction operator ends here

// Restriction starts here
Eigen::VectorXf SpaceTimeRestriction(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf Time_w_res = Time_Restriction(w_old, n, i);
    //std::cout << "\n time_w_res size: \n" << Time_w_res.size();

    Eigen::VectorXf Space_w_res = Space_Restriction(Time_w_res, i);
    std::cout << "\n w_res size: \n" << Space_w_res.size();

    return Space_w_res;
}

std::tuple<int, int, Eigen::VectorXf> Restriction(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf w_restricted = SpaceTimeRestriction(w_old, n, i);

    return std::make_tuple(n/2, i/2, w_restricted);
}

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
Eigen::VectorXf SpaceTimeProlongation(Eigen::VectorXf w_res, int n, int i){

    Eigen::VectorXf Time_w_pro = Time_Prolongation(w_res, i);
    //std::cout << "\n time_w_prolongated size: \n" << Time_w_pro.size();

    Eigen::VectorXf Space_w_pro = Space_Prolongation(Time_w_pro, i);
    std::cout << "\n w_res size: \n" << Space_w_pro.size();

    return Space_w_pro;
}

int main(){

    int n;
    int i;
    
    float f;
    
    std::cout << "Enter the number of time step: " << std::endl;

    std::cin >> n;

    std::cout << "Enter the number of space grids: " << std::endl;

    std::cin >> i;

    float delt = 1.0/n;
    float delx = 1.0/i;

    std::cout << "\n delt: " << delt << std::endl;
    std::cout << "\n delx: " << delx << std::endl;

    //Eigen::VectorXf w_old = Eigen::VectorXf::Random((i-1)*2*n);

    Eigen::VectorXf w_old = Eigen::VectorXf::Random((i-1)*2*n + 2*(i-1));

    std::cout << "\n w_old size: \n" << w_old.size() << std::endl;

    /* auto data = Restriction(w_old, n, i);
    std::cout << "\n new restricted value of n: \n" << std::get<0>(data);
    std::cout << "\n new restricted value of i: \n" << std::get<1>(data); */

    Eigen::VectorXf w_res = SpaceTimeRestriction(w_old, n, i);

    std::cout << "\n w_res size: \n" << w_res.size(); 
    std::cout << "\n new restricted value of n: \n" << n/2;
    std::cout << "\n new restricted value of i: \n" << i/2; 

    Eigen::VectorXf w_pro = SpaceTimeProlongation(w_res, (n/2), (i/2));

    
}