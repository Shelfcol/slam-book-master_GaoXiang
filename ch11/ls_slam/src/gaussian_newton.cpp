#include "gaussian_newton.h"
#include <eigen3/Eigen/Jacobi>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>

//位姿-->转换矩阵
Eigen::Matrix3d PoseToTrans(Eigen::Vector3d x)
{
    Eigen::Matrix3d trans;
    trans << cos(x(2)), -sin(x(2)), x(0),
        sin(x(2)), cos(x(2)), x(1),
        0, 0, 1;

    return trans;
}

Eigen::Matrix2d theta2R(float theta)
{
    Eigen::Matrix2d R;
    R << cos(theta), -sin(theta),
        sin(theta), cos(theta);
    return R;
}
//角度求旋转矩阵的导数
Eigen::Matrix2d theta2R_theta(float theta)
{
    Eigen::Matrix2d RT;
    RT << -sin(theta), -cos(theta),
        cos(theta), -sin(theta);
    return RT;
}

//转换矩阵－－＞位姿
Eigen::Vector3d TransToPose(Eigen::Matrix3d trans)
{
    Eigen::Vector3d pose;
    pose(0) = trans(0, 2);
    pose(1) = trans(1, 2);
    pose(2) = atan2(trans(1, 0), trans(0, 0));

    return pose;
}

//计算整个pose-graph的误差
double ComputeError(std::vector<Eigen::Vector3d> &Vertexs,
                    std::vector<Edge> &Edges)
{
    double sumError = 0;
    //对每条边进行误差求解，即里程计的值和观测值的差
    for (int i = 0; i < Edges.size(); i++)
    {
        Edge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        Eigen::Matrix3d Xi = PoseToTrans(xi);
        Eigen::Matrix3d Xj = PoseToTrans(xj);
        Eigen::Matrix3d Z = PoseToTrans(z);

        Eigen::Matrix3d Ei = Z.inverse() * Xi.inverse() * Xj;

        Eigen::Vector3d ei = TransToPose(Ei);

        sumError += ei.transpose() * infoMatrix * ei;
    }
    return sumError;
}

/**
 * @brief CalcJacobianAndError
 *         计算jacobian矩阵和error
 * @param xi    fromIdx
 * @param xj    toIdx
 * @param z     观测值:xj相对于xi的坐标
 * @param ei    计算的误差
 * @param Ai    相对于xi的Jacobian矩阵
 * @param Bi    相对于xj的Jacobian矩阵
 */
void CalcJacobianAndError(Eigen::Vector3d xi, Eigen::Vector3d xj, Eigen::Vector3d z,
                          Eigen::Vector3d &ei, Eigen::Matrix3d &Ai, Eigen::Matrix3d &Bi)
{
    //TODO--Start
    Ai.setZero();
    Bi.setZero();
    Eigen::Matrix2d Rz = theta2R(z(2));
    Eigen::Matrix2d Ri = theta2R(xi(2));
    Eigen::Matrix2d Ri_inverse = theta2R_theta(xi(2));
    Eigen::Vector2d Xij;
    Eigen::Matrix3d Z = PoseToTrans(z);

    Xij << (xj(0) - xi(0)), (xj(1) - xi(1));

    float si = sin(xi(2));
    float ci = cos(xi(2));

    // 这里的A，B是还没有乘以Rz转置的
    Ai << -ci, -si, -si * Xij(0) + ci * Xij(1),
        si, -ci, -ci * Xij(0) - si * Xij(1),
        0, 0, -1;
    Bi << ci, si, 0,
        -si, ci, 0,
        0, 0, 1;

    Eigen::Matrix3d ztinv = Z.inverse();
    ztinv(0, 2) = 0; //偏导A,B计算公式中只用了z的旋转矩阵R，所以把平移向量强制清0
    ztinv(1, 2) = 0;
    Ai = ztinv * Ai; // 乘以Z的转置，即Z的逆，这是由于旋转矩阵的关系
    Bi = ztinv * Bi;

    Eigen::Matrix3d Xi = PoseToTrans(xi);
    Eigen::Matrix3d Xj = PoseToTrans(xj);

    Eigen::Matrix3d Ei = Z.inverse() * Xi.inverse() * Xj;

    ei = TransToPose(Ei);

    //TODO--end
}

/**
 * @brief LinearizeAndSolve
 *        高斯牛顿方法的一次迭代．
 * @param Vertexs   图中的所有节点
 * @param Edges     图中的所有边
 * @return          位姿的增量
 */
Eigen::VectorXd LinearizeAndSolve(std::vector<Eigen::Vector3d> &Vertexs,
                                  std::vector<Edge> &Edges)
{
    //申请内存，H矩阵需要保存所有顶点的几个偏导值
    Eigen::MatrixXd H(Vertexs.size() * 3, Vertexs.size() * 3);
    Eigen::VectorXd b(Vertexs.size() * 3);

    H.setZero();
    b.setZero();

    //固定第一帧
    Eigen::Matrix3d I;
    I.setIdentity();
    H.block(0, 0, 3, 3) += I; //matrix.block(i,j,p,q):提取块大小为(p,q),起始于(i,j)

    //遍历每条边，构造H矩阵　＆ b向量
    for (int i = 0; i < Edges.size(); i++)
    {
        //提取信息
        Edge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        //计算误差和对应的Jacobian
        Eigen::Vector3d ei;
        Eigen::Matrix3d Ai;
        Eigen::Matrix3d Bi;

        CalcJacobianAndError(xi, xj, z, ei, Ai, Bi);

        //TODO--Start
        H.block(3 * tmpEdge.xi, 3 * tmpEdge.xi, 3, 3) += Ai.transpose() * infoMatrix * Ai;
        H.block(3 * tmpEdge.xi, 3 * tmpEdge.xj, 3, 3) += Ai.transpose() * infoMatrix * Bi;
        H.block(3 * tmpEdge.xj, 3 * tmpEdge.xi, 3, 3) += Bi.transpose() * infoMatrix * Ai;
        H.block(3 * tmpEdge.xj, 3 * tmpEdge.xj, 3, 3) += Bi.transpose() * infoMatrix * Bi;

        b.block(3 * tmpEdge.xi, 0, 3, 1) -= Ai.transpose() * infoMatrix * ei;
        b.block(3 * tmpEdge.xj, 0, 3, 1) -= Bi.transpose() * infoMatrix * ei;

        //TODO--End
    }

    //求解
    Eigen::VectorXd dx;

    //TODO--Start
    //方法1：直接求逆求解，数据量大的时候速度慢
    dx = H.inverse() * b;

    //法2：利用Eigen库求解稀疏矩阵方法
    /*
    // 默认是按列优先
    int A_rows = H.rows();
    int A_cols = H.cols();
    int b_rows = b.rows();
    int b_cols = b.cols();

    Eigen::SparseMatrix<double> A1_sparse(A_rows, A_cols);
    //Eigen::VectorXf b1_sparse(b_rows, b_cols);
   // b1_sparse = b;
    std::vector<Eigen::Triplet<double> > tripletlist;
    for (int i = 0; i < A_rows; i++)
    {
        for (int j = 0; j < A_cols; j++)
        {
            if (H(i, j) != 0)
            {
                //按Triplet方式填充，速度快
                tripletlist.push_back(Eigen::Triplet<double>(i, j, H(i, j)));

                // 直接插入速度慢
                //A1_sparse.insert(i, j) = A(i, j);
            }
        }
    }
    A1_sparse.setFromTriplets(tripletlist.begin(), tripletlist.end());

    // 压缩优化矩阵
    A1_sparse.makeCompressed();

    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > Solver_sparse;

    // 设置迭代精度
    Solver_sparse.setTolerance(0.001);
    Solver_sparse.compute(A1_sparse);

    //dx 即为解
    dx = Solver_sparse.solve(b);
    */
    //std::cout << "dx=" << dx << std::endl;


    //TODO-End

    return dx;
}
