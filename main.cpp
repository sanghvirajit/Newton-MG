#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

// Just the another way of performing Jacobi Smoother 
// Considering the whole Matrix D
Eigen::VectorXf JacobiSmoother(Eigen::MatrixXf A, Eigen::MatrixXf D, Eigen::VectorXf w_old, Eigen::VectorXf f, int NSM){

    Eigen::VectorXf w_new;
    Eigen::VectorXf w_diff;
    
    int omega=1;
    float norm;
    float tol=1E-05;

    std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << "Residual: " << "\n";
        }
        myfile.close();

    std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << "Iterations: " << "\n";
        }
        myfile_2.close();    

    for(int j=0; j < NSM; ++j){
        
        w_diff = D.inverse() * (f - A*w_old);

        w_new = w_old + omega * w_diff;

        norm = w_diff.norm();

        if(j==0){
            std::cout << "\n Residual norm: " << "\n" <<  std::endl;
        }
        
        std::cout << norm << "\n" <<  std::endl;

        std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << norm << "\n";
        }
        myfile.close();

        std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << j << "\n";
        }
        myfile_2.close();

        w_old = w_new;

    }

    std::ofstream myfile_3 ("W.txt", std::ios_base::app);

    if (myfile_3.is_open()){

            myfile_3 << w_new << "\n";
    }
    myfile_3.close();

    return w_new;

}
// Just the another way of performing Jacobi Smoother ends




// Main functions starts from here




// function to read txt file and adding values to vector elements
Eigen::VectorXf CSVtoVector(std::string file, int i, int n) {

    double data[i*n];

    Eigen::VectorXf vector = Eigen::VectorXf::Zero(i*n - 1);

    std::ifstream input(file);

    for (int j = 0; j < i*n - 1; ++j) {
        input >> data[j];
        vector[j] = data[j];

    }
    
    // std::cout << vector << std::endl;

    return vector;
}

// Saving the vector to the file
void WriteFile(Eigen::VectorXf vector, std::string filename) {
    
    std::ofstream myfile (filename.append(".txt"), std::ios_base::app);

    if (myfile.is_open()){

            myfile << vector << "\n";
    }
    myfile.close();
}

// Building Matrix A
Eigen::MatrixXf BuildA(int n,int i, float delt, float delx, float mu, float beta){

    Eigen::MatrixXf A(2*n*i + 2*i, 2*n*i + 2*i);

    //Constructing matrix D
    Eigen::MatrixXf D =  Eigen::MatrixXf::Identity(i, i);

    for(int j = 0; j < i; ++j){

        D(j, j) = -2.0;

    }

    for(int j=0; j < i-1; ++j){

        D(j+1, j) = 1.0;

    }

    for(int j=0; j < i-1; ++j){

        D(j, j+1) = 1.0;

    }

    D = 1.0/(delx*delx) * D;

    // std::cout << "\n Matrix D: \n" << D << std::endl;

    // constructing submatrix A0 for initial time step
    Eigen::MatrixXf A0 = Eigen::MatrixXf::Zero(2*i, 2*i);

    Eigen::MatrixXf a01 = Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf a02 = Eigen::MatrixXf::Zero(i, i);
    Eigen::MatrixXf a11 = -1.0 * delt * Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf a12 = (Eigen::MatrixXf::Identity(i, i)) - (delt * D);

    A0.block(0, 0, i, i) = a01;
    A0.block(0, i, i, i) = a02;
    A0.block(i, 0, i, i) = a11;
    A0.block(i, i, i, i) = a12; 

    // constructing submatrix Am for all other time steps
    Eigen::MatrixXf Am = Eigen::MatrixXf::Zero(2*i, 2*i);

    Eigen::MatrixXf Am01 = (Eigen::MatrixXf::Identity(i, i)) - (delt * D);
    Eigen::MatrixXf Am02 = (1.0/mu) * delt * Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf Am11 = -1.0 * delt * Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf Am12 = (Eigen::MatrixXf::Identity(i, i)) - (delt * D);

    Am.block(0, 0, i, i) = Am01;
    Am.block(0, i, i, i) = Am02;
    Am.block(i, 0, i, i) = Am11;
    Am.block(i, i, i, i) = Am12;

    // constructing submatrix Ant for final time step
    Eigen::MatrixXf Ant = Eigen::MatrixXf::Zero(2*i, 2*i);

    Eigen::MatrixXf Ant01 = (Eigen::MatrixXf::Identity(i, i)) - (delt * D);
    Eigen::MatrixXf Ant02 = (1.0/mu) * delt * Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf Ant11 = beta * Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf Ant12 = Eigen::MatrixXf::Identity(i, i);

    Ant.block(0, 0, i, i) = Ant01;
    Ant.block(0, i, i, i) = Ant02;
    Ant.block(i, 0, i, i) = Ant11;
    Ant.block(i, i, i, i) = Ant12;

    // constructing submatrix B
    Eigen::MatrixXf B = Eigen::MatrixXf::Zero(2*i, 2*i);
    Eigen::MatrixXf B01 = (-1.0) * Eigen::MatrixXf::Identity(i, i);
    B.block(0, 0, i, i) = B01;

    // constructing submatrix C
    Eigen::MatrixXf C = Eigen::MatrixXf::Zero(2*i, 2*i);
    Eigen::MatrixXf C11 = (-1.0) * Eigen::MatrixXf::Identity(i, i);
    C.block(i, i, i, i) = C11;

    // substituting all the submatrices to matrix A
    A.block(0, 0, 2*i, 2*i) = A0;
    A.block(2*n*i, 2*n*i, 2*i, 2*i) = Ant;

    for(int j=1; j < n ; ++j){

        A.block(2*i*j, 2*i*j, 2*i, 2*i) = Am;
        
    } 

    for(int j=0; j < n-1; ++j){

        A.block(2*i*(1+j), 2*i*j, 2*i, 2*i) = B;
    }
        

    for(int j=0; j < n-1; ++j){

        A.block(2*i*j, 2*i*(1+j), 2*i, 2*i) = C;
    } 

    return A;
} 
// Building Matrix A ends


// Building vector f
Eigen::VectorXf BoundaryConditions(Eigen::VectorXf f, Eigen::VectorXf Y, Eigen::VectorXf Z, Eigen::VectorXf P, int i, int n, float delt, float beta){

    f.block(0, 0, i, 1) = Y.block(0, 0, i, 1);
        
    for(int j=0; j < n; ++j){

        f.block(2*i*j + i, 0, i, 1) = -1.0 * delt * Z.block(i*j, 0, i, 1);

    }
    
    f.block(2*i*n + i, 0, i, 1) = -1.0 * beta * delt * Z.block(i*n, 0, i, 1);

    return f;
}
// Building vector f ends


// Performing LDU decomposition
std::vector<Eigen::MatrixXf> LUdecomposition(Eigen::MatrixXf M){

    std::vector<Eigen::MatrixXf> myVec;

    Eigen::MatrixXf U = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    U.triangularView<Eigen::StrictlyUpper>() = M.triangularView<Eigen::StrictlyUpper>();
    
    Eigen::MatrixXf L = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    L.triangularView<Eigen::StrictlyLower>() = M.triangularView<Eigen::StrictlyLower>();
    
    Eigen::MatrixXf D = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    for(int i=0; i < M.rows(); ++i){

        D(i, i) = M(i, i);

    }
    
    myVec.push_back(L);
    myVec.push_back(D);
    myVec.push_back(U);

    return myVec;
}
// Performing LDU decomposition ends


// Building Jacobi Smoother
// Solving in sequential manner for each time step
Eigen::VectorXf SpaceTimeBlockJacobiSmoother(Eigen::MatrixXf A, Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, int Iterations, float Tol){

    Eigen::VectorXf d;
    Eigen::VectorXf w_new;
    Eigen::VectorXf w_diff;

    float norm = 1.0;
    float omega = 0.01;

    std::vector<Eigen::MatrixXf> myVec= LUdecomposition(A);
    Eigen::MatrixXf L = myVec[0];
    Eigen::MatrixXf D = myVec[1];
    Eigen::MatrixXf U = myVec[2];

    // adding initial line to the file
    std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << "Residual: " << "\n";
        }
        myfile.close();

    // adding initial line to the file
    std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << "Iterations: " << "\n";
        }
        myfile_2.close();    

    int j = 0;
    while( (norm > Tol) && (j < Iterations) ){
        
        // calculating the defect
        d = f - A * w_old;

        for(int k=0; k < n; ++k){
            
            d.block(2*i*k, 0, 2*i, 1) = D.block(2*i*k, 2*i*k, 2*i, 2*i).inverse() * d.block(2*i*k, 0, 2*i, 1);
              
        }

        // Update the vector
        w_new = w_old + omega * d;

        // calculating the norm
        w_diff = w_new - w_old;

        //std::cout << "\n w_diff: \n" << w_diff;

        norm = w_diff.norm();

        // Printing the residual norm
        if(j==0){
            std::cout << "\n Residual norm: " << "\n" <<  std::endl;
        }
        std::cout << norm << "\n" <<  std::endl;
        
        // writing residual norm to the file
        std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << norm << "\n";
        }
        myfile.close();

        // writing number of iterations to the file
        std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << j << "\n";
        }
        myfile_2.close();

        w_old = w_new;
        ++j;
    }
    std::cout << "\n Solution Converged in: "  << j << " Iterations." << std::endl;

    return w_new;
}
// Building Jacobi Smoother ends


// Building Gauss Seidel Smoother
Eigen::VectorXf SpaceTimeBlockGaussSeidelSmoother(Eigen::MatrixXf A, Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, int Iterations, float Tol){

    Eigen::VectorXf d;
    Eigen::VectorXf w_new;
    Eigen::VectorXf w_diff;

    float norm = 1;
    float omega = 0.01;

    std::vector<Eigen::MatrixXf> myVec= LUdecomposition(A);
    Eigen::MatrixXf L = myVec[0];
    Eigen::MatrixXf D = myVec[1];
    Eigen::MatrixXf U = myVec[2];

    // adding initial line to the file
    std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << "Residual: " << "\n";
        }
        myfile.close();

    // adding initial line to the file
    std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << "Iterations: " << "\n";
        }
        myfile_2.close();    

    int j = 0;
    while( (j < Iterations) && (norm > Tol) ){
        
        // calculating the defect
        d = f - A * w_old;

        for(int k=0; k < n; ++k){
            
            d.block(2*i*k, 0, 2*i, 1) = ( L.block(2*i*k, 2*i*k, 2*i, 2*i) + D.block(2*i*k, 2*i*k, 2*i, 2*i)).inverse() * d.block(2*i*k, 0, 2*i, 1);
              
        }

        w_new = w_old + omega * d;
        
        // calculating the norm
        w_diff = w_new - w_old;

        norm = w_diff.norm();

        // Printing the norm
        if(j==0){
            std::cout << "\n Residual norm: " << "\n" <<  std::endl;
        }
        std::cout << norm << "\n" <<  std::endl;

        // Writing the residual norm to the file 
        std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << norm << "\n";
        }
        myfile.close();

        // Writing number of iterations to the file 
        std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << j << "\n";
        }
        myfile_2.close();

        w_old = w_new;

        ++j;
    }
    
    std::cout << "\n Solution Converged in: "  << j << " Iterations " << std::endl;

    return w_new;

}
// Building Gauss Seidel Smoother ends

// Building SOR Smoother
Eigen::VectorXf SpaceTimeBlockSORSmoother(Eigen::MatrixXf A, Eigen::VectorXf w_old, Eigen::VectorXf f, float omega_2, int n, int i, int Iterations, float Tol){

    Eigen::VectorXf d;
    Eigen::VectorXf w_new;
    Eigen::VectorXf w_diff;

    float norm = 1;
    float omega_1 = 0.03;

    std::vector<Eigen::MatrixXf> myVec= LUdecomposition(A);
    Eigen::MatrixXf L = myVec[0];
    Eigen::MatrixXf D = myVec[1];
    Eigen::MatrixXf U = myVec[2];

    // adding initial line to the file
    std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << "Residual: " << "\n";
        }
        myfile.close();

    // adding initial line to the file
    std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << "Iterations: " << "\n";
        }
        myfile_2.close();    

    int j = 0;
    while( (j < Iterations) && (norm > Tol) ){
        
        // calculating the defect
        d = f - A * w_old;

        for(int k=0; k < n; ++k){
            
            d.block(2*i*k, 0, 2*i, 1) = ( L.block(2*i*k, 2*i*k, 2*i, 2*i) + (1/omega_2) * D.block(2*i*k, 2*i*k, 2*i, 2*i)).inverse() * d.block(2*i*k, 0, 2*i, 1);
              
        }

        w_new = w_old + omega_1 * d;
        
        // calculating the norm
        w_diff = w_new - w_old;

        norm = w_diff.norm();

        // Printing the norm
        if(j==0){
            std::cout << "\n Residual norm: " << "\n" <<  std::endl;
        }
        std::cout << norm << "\n" <<  std::endl;

        // Writing the residual norm to the file 
        std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << norm << "\n";
        }
        myfile.close();

        // Writing number of iterations to the file 
        std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << j << "\n";
        }
        myfile_2.close();

        w_old = w_new;

        ++j;
    }
    
    std::cout << "\n Solution Converged in: "  << j << " Iterations " << std::endl;

    return w_new;

}
// Building SOR Smoother ends

// Restriction starts here

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

    //std::cout << "\n restriction operator: \n" << time_block;

    for(int j=0; j < n; ++j){

        space_res_operator.block((i/2 - 1) * 2 * j, (i-1) * 2 * j, (i/2 - 1) * 2, (i-1) * 2) = time_block;
    }

    //std::cout << "\n space_res_operator: \n" << space_res_operator;

    d_space_restricted = (1.0/4) * space_res_operator * d;

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
    std::cout << "\n w_res size: " << Space_w_res.size();

    return Space_w_res;
}

std::tuple<int, int, Eigen::VectorXf> Restriction(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf w_restricted = SpaceTimeRestriction(w_old, n, i);

    return std::make_tuple(n/2, i/2, w_restricted);
}
// Restriction ends here


// Prolongation starts here
// Space Prolongation starts here
Eigen::VectorXf Space_Prolongation(Eigen::VectorXf d, int i){

    Eigen::VectorXf d_space_prolongated;
    
    int n = d.size()/((i-1)*2);

    Eigen::MatrixXf space_prolongation = Eigen::MatrixXf::Zero( (2*(i-1) + 1) * 2 * n, d.size());

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

    //std::cout << "\n Prolongation_operator: \n" << Prolongation_operator;

    for(int j=0; j < n; ++j){

        space_prolongation.block((2*(i-1) + 1)* 2 * j, (i-1) * 2 * j, (2*(i-1) + 1)* 2, (i-1) * 2) = Prolongation_operator;
    }

    // std::cout << "\n space_prolongation: \n" << space_prolongation;

    d_space_prolongated = (1.0/4) * space_prolongation * d;

    return d_space_prolongated;
}
//Space Prolongation ends here

// Time prolongation starts here
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
// Time prolongation ends here

// Combining time and space prolongation
Eigen::VectorXf SpaceTimeProlongation(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf Time_w_pro = Time_Prolongation(w_old, i);
    //std::cout << "\n time_w_prolongated size: \n" << Time_w_pro.size();

    Eigen::VectorXf Space_w_pro = Space_Prolongation(Time_w_pro, i);
    std::cout << "\n w_prol size: " << Space_w_pro.size();

    return Space_w_pro;
}

std::tuple<int, int, Eigen::VectorXf> Prolongation(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf w_prol = SpaceTimeProlongation(w_old, n, i);

    return std::make_tuple(n*2, i*2, w_prol);
}
// Prolongation ends here

// Single cycle starts here
Eigen::VectorXf CoarseGridSolver(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, float beta){
    
    float delt = 1.0/n;
    float delx = 1.0/i;

    Eigen::MatrixXf A = BuildA(n, (i-1), delt, delx, mu, beta);
    
    const clock_t begin_time = clock();
    Eigen::VectorXf w_tilda = SpaceTimeBlockJacobiSmoother(A, w_old, f, n, (i-1), 20000, 1E-08);
    std::cout << "\n Execution time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC;

    return w_tilda;

}
// Single cycle ends here

// Single cycle starts here
Eigen::VectorXf SingleCycle(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, float beta){
    
    float delt = 1.0/n;
    float delx = 1.0/i;

    Eigen::MatrixXf A = BuildA(n, (i-1), delt, delx, mu, beta);
    
    const clock_t begin_time = clock();
    Eigen::VectorXf w = SpaceTimeBlockSORSmoother(A, w_old, f, 1.5, n, (i-1), 20000, 2.85E-02);
    std::cout << "\n Execution time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC;

    WriteFile(w, "w_final");

    return w;

}
// Single cycle ends here

// Two Grid V-cycle starts here
Eigen::VectorXf VCycle(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, float beta){
    
    remove("w_tilda.txt");
    remove("w_tilda_vector.txt");
    remove("w_final.txt");

    float delt = 1.0/n;
    float delx = 1.0/i;
    
    Eigen::MatrixXf A = BuildA(n, (i-1), delt, delx, mu, beta);
    
    Eigen::VectorXf w_tilda = SpaceTimeBlockJacobiSmoother(A, w_old, f, n, (i-1), 500, 1E-04);

    WriteFile(w_tilda, "w_tilda");
    
    Eigen::VectorXf d = f - A * w_tilda;

    std::cout << "\n \n Restricting.... \n" << std::endl;

    auto restriction_data = Restriction(d, n, i);
    n = std::get<0>(restriction_data);
    i = std::get<1>(restriction_data);
    delt = 1.0/n;
    delx = 1.0/i;

    std::cout << "\n Restricted value of n: " << n << std::endl;
    std::cout << "\n Restricted value of i: " << i << std::endl;

    std::cout << "\n delt: " << delt << std::endl;
    std::cout << "\n delx: " << delx << std::endl;

    Eigen::MatrixXf A_2h = BuildA(n, (i-1), delt, delx, mu, beta);

    Eigen::VectorXf error = Eigen::VectorXf::Zero(std::get<2>(restriction_data).size());
    
    std::cout << "\n Solving residual equation on course grid..... \n" << std::endl;

    // Course grid solver
    Eigen::VectorXf d_tilda_tilda = CoarseGridSolver(error, std::get<2>(restriction_data), n, i, mu, beta);
    
    std::cout << "\n \n Prolongating.... \n" << std::endl;

    auto prolongation_data = Prolongation(d_tilda_tilda, n, i);
    n = std::get<0>(prolongation_data);
    i = std::get<1>(prolongation_data);
    delt = 1.0/n;
    delx = 1.0/i;
    
    std::cout << "\n Prolongated value of n: " << n << std::endl;
    std::cout << "\n Prolongated value of i: " << i << std::endl;

    std::cout << "\n delt: " << delt << std::endl;
    std::cout << "\n delx: " << delx << std::endl; 

    Eigen::VectorXf w_tilda_vector = w_tilda + std::get<2>(prolongation_data); 

    WriteFile(w_tilda_vector, "w_tilda_vector");

    //Post smoothing
    Eigen::VectorXf w = SpaceTimeBlockJacobiSmoother(A, w_tilda_vector, f, n, (i-1), 500, 1E-04); 

    WriteFile(w, "w_final");

    return w;
}
// Two Grid V-cycle starts here

// Multilevel V cycle starts here
Eigen::VectorXf SpaceTimeMultigrid(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, float beta, int K){

    Eigen::MatrixXf A;
    Eigen::VectorXf w_tilda, w_tilda_vector, w;
    Eigen::VectorXf error, d, residual;

    float delt = 1.0/n;
    float delx = 1.0/i;
    
    A = BuildA(n, (i-1), delt, delx, mu, beta);

    if (K==1){

        return SpaceTimeBlockSORSmoother(A, w_old, f, 1.5, n, (i-1), 20000, 1E-05);
    
    }
    
    w_tilda = SpaceTimeBlockSORSmoother(A, w_old, f, 1.5, n, (i-1), 10, 1E-08);   
   
    d = f - A * w_tilda;

    std::cout << "\n \n Restricting.... \n" << std::endl;

    auto d_tilda = Restriction(d, n, i);

    n = std::get<0>(d_tilda);
    i = std::get<1>(d_tilda);
    delt = 1.0/n;
    delx = 1.0/i;

    Eigen::MatrixXf A_2h = BuildA(n, (i-1), delt, delx, mu, beta);

    error = SpaceTimeMultigrid(Eigen::VectorXf::Zero(std::get<2>(d_tilda).size()), std::get<2>(d_tilda), n, i, mu, beta, K-1);

    std::cout << "\n \n Prolongating.... \n" << std::endl;

    auto d_tilda_tilda_vector = Prolongation(error, std::get<0>(d_tilda), std::get<1>(d_tilda));

    w_tilda_vector =  w_tilda + std::get<2>(d_tilda_tilda_vector);
    
    w = SpaceTimeBlockSORSmoother(A, w_old, f, 1.5, n, (i-1), 10, 1E-08);   
    
    return w;
}
// // Multilevel V cycle ends here

// Multiple Multilevel V cycle starts here
Eigen::VectorXf Multiple_SpaceTimeMultigrid(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, float beta, int K){

    Eigen::MatrixXf A;
    Eigen::VectorXf w_tilda, w_tilda_vector, w;
    Eigen::VectorXf error, d, residual;

    float residual_norm = 1;
    float residual_norm_2 = 1;
    float tol = 0.001;

    float delt = 1.0/n;
    float delx = 1.0/i;
    int j = 0;

    A = BuildA(n, (i-1), delt, delx, mu, beta);

    if (K==1){

        return SpaceTimeBlockSORSmoother(A, w_old, f, 1.5, n, (i-1), 20000, 1E-05);

    }
    
    while((residual_norm >= tol) && ((residual_norm_2 >= tol))){

        std::cout << "\n Multigrid V cycle Converged in: "  << j+1 << " Cycles " << std::endl;
        std::cout << "\n Residual norm (f - Aw): "  << residual_norm << std::endl;
        std::cout << "\n Residual norm (w - w_old): "  << residual_norm_2 << std::endl;

        w_tilda = SpaceTimeBlockSORSmoother(A, w_old, f, 1.5, n, (i-1), 15, 1E-08);   

        d = f - A * w_tilda;

        std::cout << "\n \n Restricting.... \n" << std::endl;

        auto d_tilda = Restriction(d, n, i);

        float R_n = std::get<0>(d_tilda);
        float R_i = std::get<1>(d_tilda);
        delt = 1.0/R_n;
        delx = 1.0/R_i;

        Eigen::MatrixXf A_2h = BuildA(R_n, (R_i-1), delt, delx, mu, beta);

        error = SpaceTimeMultigrid(Eigen::VectorXf::Zero(std::get<2>(d_tilda).size()), std::get<2>(d_tilda), R_n, R_i, mu, beta, K-1);

        std::cout << "\n \n Prolongating.... \n" << std::endl;

        auto d_tilda_tilda_vector = Prolongation(error, std::get<0>(d_tilda), std::get<1>(d_tilda));
        
        w_tilda_vector =  w_tilda + std::get<2>(d_tilda_tilda_vector);

        w = SpaceTimeBlockSORSmoother(A, w_tilda_vector, f, 1.5, n, (i-1), 15, 1E-08);

        residual = f - A * w;

        residual_norm = residual.norm();

        ++j;

        residual_norm_2 = (w - w_old).norm();

        w_old = w;

    }

    std::cout << "\n Multigrid V cycle Converged in: "  << j+1 << " Cycles " << std::endl;
    std::cout << "\n Residual norm (f - Aw): "  << residual_norm << std::endl;
    std::cout << "\n Residual norm (w - w_old): "  << residual_norm_2 << std::endl;

    WriteFile(w_tilda, "w_tilda");
    WriteFile(w_tilda_vector, "w_tilda_vector");
    WriteFile(w, "w_final");

    return w;
}
// Multiple Multilevel V cycle ends here

int main(){

    int n;
    int i;
    
    float mu = 1;
    float beta = 0;

    remove("Residual.txt");
    remove("Iterations.txt");
    remove("w_tilda.txt");
    remove("w_tilda_vector.txt");
    remove("w_final.txt");

    std::cout << "Enter the number of time step: " << std::endl;

    std::cin >> n;

    std::cout << "Enter the number of space grids: " << std::endl;

    std::cin >> i;

    float delt = 1.0/n;
    float delx = 1.0/i;

    std::cout << "\n delt: " << delt << std::endl;
    std::cout << "\n delx: " << delx << std::endl;

    Eigen::VectorXf Y = CSVtoVector("Y.txt", i, n);
    Eigen::VectorXf P = CSVtoVector("P.txt", i, n);
    Eigen::VectorXf Z = CSVtoVector("Z.txt", i, n); 

    Eigen::VectorXf w_old = Eigen::VectorXf::Random(2*(i-1)*n + 2*(i-1));
    
    Eigen::VectorXf f = Eigen::VectorXf::Zero(2*(i-1)*n + 2*(i-1)); 
    f = BoundaryConditions(f, Y, Z, P, i-1, n, delt, beta);
    // std::cout << " \n Vector f: \n" << f << std::endl;
    
    // const clock_t begin_time = clock();
    // Eigen::VectorXf w = SingleCycle(w_old, f, n, i, mu, beta);
    // std::cout << "\n Execution time: " << float( clock () - begin_time ) / CLOCKS_PER_SEC;

    const clock_t begin_time = clock();
    Eigen::VectorXf w = Multiple_SpaceTimeMultigrid(w_old, f, n, i, mu, beta, 3);
    std::cout << "\n Execution time: " << float( clock () - begin_time ) / CLOCKS_PER_SEC;

    return 0;
}  