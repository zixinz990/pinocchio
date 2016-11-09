#ifndef __se3_ldlt6_hpp__
#define __se3_ldlt6_hpp__

#include "pinocchio/spatial/inertia.hpp"


namespace se3
{

  class LDLT6
  {
  public:
    typedef Inertia::Symmetric3 Symmetric3;
    typedef Inertia::Matrix6 Matrix6;
    typedef Inertia::Vector6 Vector6;

    typedef const Eigen::TriangularView<const Matrix6, Eigen::UnitLower> MatrixL;
    typedef const Eigen::Transpose<const Matrix6> Matrix6T;
    typedef const Eigen::TriangularView<const Matrix6T, Eigen::UnitUpper> MatrixU;  
    typedef const Eigen::DiagonalWrapper<const Vector6> Diag6;

  public: // --- Builders

    // WIP constructor. Should be set to NAN in release.
    LDLT6() : L(Matrix6::Identity()),D(Vector6::Ones()) {}

    /*
      d0  = y0
      l10 = y1/y0
      l20 = y2/y0
      d1  = (y3 - l10*y1)
      l21 = (y4-l20*y1)/d1
      d2  = (y5 - l20*y2 - d1*l21**2)
    */

    template<typename DM, typename DV>
    static void chol3( const Symmetric3& S,
                       Eigen::MatrixBase<DM> & L,
                       Eigen::MatrixBase<DV> & D)
    {
      const Vector6 & y = S.data();
      D(0,0) =  y[0]         ;
      L(1,0) =  y[1] / y[0]  ;
      L(2,0) =  y[3] / y[0]  ;
      D(1,0) =  y[2] - L(1,0) * y[1] ;
      L(2,1) = (y[4] - L(2,0) * y[1] ) / D[1];
      D(2,0) =  y[5] - L(2,0) * y[3] - D[1]*L(2,1)*L(2,1);
    }

    void fromInertia( const Inertia & Y )
    {
      L(1,0) = L(2,0) = L(2,1) = 0.0;
      D.head<3>().fill(Y.mass());
    
      L.bottomLeftCorner<3,3>() = skew( Y.lever() );

      Eigen::Block<Matrix6,3,3> Lr = L.bottomRightCorner<3,3>();
      Eigen::VectorBlock<Vector6,3> Dr = D.tail<3>();
      chol3( Y.inertia(), Lr,Dr );
    }


  public: // --- Utils


  public: // --- Debug

    MatrixL matrixL() const { return L; }
    MatrixU matrixU() const { return (Matrix6T)L.transpose(); }
    Diag6   matrixD() const { return D; }

    // Matrix6 matrixL() const
    // {
    //   Matrix6 Ldense;
    //   Ldense = Matrix6::Identity();
    //   for(int r=0;r<6;++r )
    //     for(int c=0;c<r;++c )
    //       Ldense(r,c) = L(r,c);

    //   return Ldense;
    // }
    // Matrix6 matrixD() const
    // {
    //   return (Matrix6)D.asDiagonal();
    // }

    // Reconstruct the matrix L*D*L.T
    Matrix6 reconstruct()
    {
      return matrixL()*Matrix6(matrixD())*matrixU();
    }

  public: // protected:
    Matrix6 L;
    Vector6 D;

  };

} // namespace se3

#endif // #ifndef __se3_ldlt6_hpp__
