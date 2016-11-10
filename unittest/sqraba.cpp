//
// Copyright (c) 2016 CNRS
//
// This file is part of Pinocchio
// Pinocchio is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// Pinocchio is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// Pinocchio If not, see
// <http://www.gnu.org/licenses/>.

#include "pinocchio/spatial/fwd.hpp"
#include "pinocchio/spatial/se3.hpp"
#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/sqraba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/spatial/ldlt6.hpp"

#include "pinocchio/algorithm/compute-all-terms.hpp"
#include "pinocchio/tools/timer.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>
#include <Eigen/Cholesky>

BOOST_AUTO_TEST_SUITE ( BOOST_TEST_MODULE )


void L_times_x( const se3::Inertia::Matrix6& L,
               const se3::Inertia::Vector6& x,
               se3::Inertia::Vector6& Lx )
{
  //for i in range(6):
  //print "Lx[%d] =" % i, "" if i==0 else " L(%d,%d)*x[%d] +"*(i) % reduce(lambda x,y:x+y,[ (i,j,j) for j in range(i) ]), "x[%d];" % (i)

  Lx[0] =  x[0];
  Lx[1] =  L(1,0)*x[0] + x[1];
  Lx[2] =  L(2,0)*x[0] + L(2,1)*x[1] + x[2];
  Lx[3] =  L(3,0)*x[0] + L(3,1)*x[1] + L(3,2)*x[2] + x[3];
  Lx[4] =  L(4,0)*x[0] + L(4,1)*x[1] + L(4,2)*x[2] + L(4,3)*x[3] + x[4];
  Lx[5] =  L(5,0)*x[0] + L(5,1)*x[1] + L(5,2)*x[2] + L(5,3)*x[3] + L(5,4)*x[4] + x[5];
}

void U_times_x( const se3::Inertia::Matrix6& L,
                const se3::Inertia::Vector6& x,
                se3::Inertia::Vector6& Ux )
{
  // for i in range(6):
  //     print             \
  //         "Ux[%d] =" % i,                    \
  //         "x[%d]" % i,                           \
  //         "" if i==5 else "+ L(%d,%d)*x[%d]"*(5-i)                   \
  //         % reduce(lambda x,y:x+y,[ (i+j+1,i,i+j+1) for j in range(5-i) ]), \
  //         ";"

  Ux[0] = x[0] + L(1,0)*x[1]+ L(2,0)*x[2]+ L(3,0)*x[3]+ L(4,0)*x[4]+ L(5,0)*x[5] ;
  Ux[1] = x[1] + L(2,1)*x[2]+ L(3,1)*x[3]+ L(4,1)*x[4]+ L(5,1)*x[5] ;
  Ux[2] = x[2] + L(3,2)*x[3]+ L(4,2)*x[4]+ L(5,2)*x[5] ;
  Ux[3] = x[3] + L(4,3)*x[4]+ L(5,3)*x[5] ;
  Ux[4] = x[4] + L(5,4)*x[5] ;
  Ux[5] = x[5]  ;
}

void D_times_x( const se3::Inertia::Vector6& D,
                const se3::Inertia::Vector6& x,
                se3::Inertia::Vector6& Lx )
{
  //for i in range(6):
  //print "Lx[%d] = D[%d] * x[%d];" % (i,i,i)
  Lx[0] = D[0] * x[0];
  Lx[1] = D[1] * x[1];
  Lx[2] = D[2] * x[2];
  Lx[3] = D[3] * x[3];
  Lx[4] = D[4] * x[4];
  Lx[5] = D[5] * x[5];
}


               

BOOST_AUTO_TEST_CASE ( test_ldlt )
{
  typedef se3::Inertia::Matrix6 Matrix6;
  typedef se3::Inertia::Vector6 Vector6;

  se3::Inertia Y = se3::Inertia::Random();
  se3::LDLT6 ldlt; 
  ldlt.fromInertia(Y);

  BOOST_CHECK( Y.matrix().isApprox( ldlt.reconstruct())  );

  Eigen::LDLT<se3::Inertia::Matrix6> ldltX(Y.matrix());

  // typedef se3::LDLT6::MatrixL MatrixL;
  // typedef se3::LDLT6::MatrixU MatrixU;
  // typedef se3::LDLT6::Matrix6T Matrix6T;
  // typedef se3::LDLT6::Diag6 Diag6;

  std::cout << "L = " << (Matrix6)ldlt.matrixL() << std::endl<< std::endl<< std::endl;
  std::cout << "D = " << (Matrix6)ldlt.matrixD() << std::endl<< std::endl<< std::endl;
  std::cout << "U = " << (Matrix6)ldlt.matrixU() << std::endl<< std::endl<< std::endl;
  
  Vector6 x = Vector6::Random();
  Vector6 y;

  L_times_x(ldlt.L,x,y);
  BOOST_CHECK( y.isApprox( ldlt.matrixL()*x) );
  BOOST_CHECK( (ldlt.L*x).isApprox( ldlt.matrixL()*x) );

  U_times_x(ldlt.L,x,y);
  BOOST_CHECK( y.isApprox( ldlt.matrixU()*x) );
  BOOST_CHECK( (ldlt.L.transpose()*x).isApprox( ldlt.matrixU()*x) );

  D_times_x(ldlt.D,x,y);
  BOOST_CHECK( y.isApprox( ldlt.matrixD()*x) );
}

BOOST_AUTO_TEST_CASE( test_perf_ldlt )
{
  StackTicToc timer(StackTicToc::US); 
  const int NBT = 1000000;

  typedef se3::Inertia::Matrix6 Matrix6;
  typedef se3::Inertia::Vector6 Vector6;
  typedef Eigen::LDLT<Matrix6> LDLTx6;

  se3::container::aligned_vector<se3::Inertia> Ys(NBT);
  se3::container::aligned_vector<se3::LDLT6  > sYs(NBT);
  se3::container::aligned_vector<     Matrix6> Ms(NBT);
  se3::container::aligned_vector<     LDLTx6 > sMs(NBT);  
  se3::container::aligned_vector<     Vector6> xs(NBT);  
  se3::container::aligned_vector<     Vector6> ys(NBT);  
  
  for(int i=0;i<NBT;++i)
    {
      Ys[i] = se3::Inertia::Random();
      Ms[i] = Ys[i].matrix();
      xs[i] = Vector6::Random();
      ys[i] = Vector6::Random();
    }

  double force  = 0.0;
  int    sample = int(Ms[55](0,0)/10*NBT);
  std::cout << "*** SAMPLE = " << sample << std::endl;

  
  std::cout << "Measuring tailored LDLT = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      sYs[_smooth].fromInertia(Ys[_smooth]);
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  sYs[sample].L(1,0) << std::endl;

  std::cout << "Measuring Eigen LDLT = \t\t";
  timer.tic();
  SMOOTH(NBT)
    {
      sMs[_smooth].compute(Ms[_smooth]);
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  sMs[sample].matrixL()(1,0) << std::endl;

  //std::cout << std::endl << std::endl << Matrix6(sYs[40].matrixL()) << std::endl;
  //std::cout << std::endl << std::endl << Matrix6(sMs[40].matrixL()) << std::endl;

  std::cout << std::endl;
  std::cout << "Measuring Dense L product = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      ys[_smooth] = sYs[_smooth].L * xs[_smooth];
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  ys[sample][1] << std::endl;

  std::cout << "Measuring sparse L1 product = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      ys[_smooth] = sYs[_smooth].matrixL() * xs[_smooth];
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  ys[sample][1] << std::endl;

  std::cout << "Measuring tailored L product = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      L_times_x( sYs[_smooth].L, xs[_smooth], ys[_smooth]);
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  ys[sample][1] << std::endl;

  std::cout << std::endl;
  std::cout << "Measuring Dense U product = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      ys[_smooth] = sYs[_smooth].L.transpose() * xs[_smooth];
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  ys[sample][1] << std::endl;

  std::cout << "Measuring sparse U1 product = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      ys[_smooth] = sYs[_smooth].matrixU() * xs[_smooth];
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  ys[sample][1] << std::endl;

  std::cout << "Measuring tailored U1 product = ";
  timer.tic();
  SMOOTH(NBT)
    {
      U_times_x(sYs[_smooth].L,xs[_smooth],ys[_smooth]);
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  ys[sample][1] << std::endl;

  std::cout << std::endl;
  std::cout << "Measuring Sparse D product = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      ys[_smooth] = sYs[_smooth].matrixD() * xs[_smooth];
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  ys[sample][1] << std::endl;

  std::cout << "Measuring tailored D product = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      D_times_x(sYs[_smooth].D,xs[_smooth],ys[_smooth]);
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  ys[sample][1] << std::endl;

  /* --- CHOL 3 --- */
  typedef Eigen::Matrix3d Matrix3;
  typedef Eigen::Vector3d Vector3;
  typedef Eigen::LDLT<Matrix3> LDLT3;

  se3::container::aligned_vector<se3::Symmetric3> S3(NBT);  
  se3::container::aligned_vector<     Matrix3   > L3(NBT);  
  se3::container::aligned_vector<     Vector3   > D3(NBT);  
  se3::container::aligned_vector<     LDLT3     > ldlt3(NBT);  
  se3::container::aligned_vector<     Matrix3   > M3(NBT);  
  for(int i=0;i<NBT;++i)
    {
      S3[i] = se3::Symmetric3::Random();
      M3[i] = S3[i].matrix();
    }

  std::cout << std::endl;
  std::cout << "Measuring tailored chol 3 = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      se3::LDLT6::chol3(S3[_smooth],L3[_smooth],D3[_smooth]);
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  D3[sample][1] << std::endl;

  std::cout << "Measuring eigen chol 3 = \t";
  timer.tic();
  SMOOTH(NBT)
    {
      ldlt3[_smooth].compute( M3[_smooth]);
    }
  timer.toc(std::cout,NBT);
  std::cout << "\t\t\t\t\t\t\tForce check" <<  ldlt3[sample].matrixL()(1,0) << std::endl;

}

BOOST_AUTO_TEST_CASE ( test_sqraba_simple )
{
  using namespace Eigen;
  using namespace se3;

  se3::Model model; buildModels::humanoidSimple(model);
  
  se3::Data data(model);
  se3::Data data_ref(model);

  VectorXd q = VectorXd::Ones(model.nq);
  VectorXd v = VectorXd::Ones(model.nv);
  VectorXd tau = VectorXd::Zero(model.nv);
  VectorXd a = VectorXd::Ones(model.nv);
  
  computeAllTerms(model, data_ref, q, v);
 
  se3::Inertia & Y = model.inertias[1];
  se3::LDLT6 ldlt; 
  ldlt.fromInertia(Y);
  BOOST_CHECK( Y.matrix().isApprox( ldlt.reconstruct())  );

  // std::cout << "Y = " << Y.matrix() << std::endl;
  // std::cout << "L = " << ldlt.matrixL() << std::endl;
  // std::cout << "D = " << ldlt.matrixD() << std::endl;
  
}

BOOST_AUTO_TEST_SUITE_END ()
