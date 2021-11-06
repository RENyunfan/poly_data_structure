//
// Created by yunfan on 2021/11/6.
//

#ifndef POLY_DATA_STRUCTURE_PIECE_HPP
#define POLY_DATA_STRUCTURE_PIECE_HPP
#include "common_include.hpp"
// A single piece of a trajectory, which is indeed a polynomial
class Piece {
private:
    // Piece(t) = c5*t^5 + c4*t^4 + ... + c1*t + c0
    // The natural coefficient matrix = [c5,c4,c3,c2,c1,c0]
    double duration{0};
    // Any time in [0, T] is normalized into [0.0, 1.0]
    // Therefore, nCoeffMat = [c5*T^5,c4*T^4,c3*T^3,c2*T^2,c1*T,c0*1]
    // is used for better numerical stability
    //    CoefficientMat nCoeffMat;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nCoeffMat;
    int order_;
    bool had_sampled_ = false;
    std::vector<Vec3> sampled_positions_;

    bool had_length_ = false;
    double length_;
    bool is_empty_{true};

    double cost_;

    DynamicMat fw_mat_;
    bool has_fw_mat_ = false;

public:
    Piece() = default;

    // Constructor from duration and coefficient
    Piece(double dur, DynamicMat coeffs) : duration(dur) {
        order_ = coeffs.cols() - 1;
        nCoeffMat.resize(3, coeffs.cols());
        double t = 1.0;
        for (int i = order_; i >= 0; i--) {
            nCoeffMat.col(i) = coeffs.col(i) * t;
            t *= dur;
        }
        is_empty_ = false;
    }

    // Constructor from boundary condition and duration
    Piece(DynamicMat boundCond, double dur) : duration(dur) {

        if (boundCond.cols() == 8) {
            order_ = 7;
            nCoeffMat.resize(3, order_ + 1);
            // The BoundaryCond matrix boundCond = [p(0),v(0),a(0),p(T),v(T),a(T)]
            /*                                          0   1    2    3    4   5     6    7
             *   The BoundaryCond matrix boundCond = [p(0),v(0),a(0),j(0),p(T),v(T),a(T),j(T)]
             * */
            double t1 = dur;
            double t2 = t1 * t1;
            double t3 = t2 * t1;


            // Inverse mapping is computed without explicit matrix inverse
            // It maps boundary condition to normalized coefficient matrix
            Eigen::Array3d p_0 = boundCond.col(0),
            v_0 = boundCond.col(1),
            a_0 = boundCond.col(2),
            j_0 = boundCond.col(3),
            p_f = boundCond.col(4),
            v_f = boundCond.col(5),
            a_f = boundCond.col(6),
            j_f = boundCond.col(7);

            nCoeffMat.col(0) =
                    20 * p_0 - 20 * p_f + 10 * t1 * v_0 + 10 * t1 * v_f + 2 * t2 * a_0 - 2 * t2 * a_f + (t3 * j_0) / 6 +
                    (t3 * j_f) / 6;
            nCoeffMat.col(1) =
                    70 * p_f - 70 * p_0 - 36 * t1 * v_0 - 34 * t1 * v_f - (15 * t2 * a_0) / 2 + (13 * t2 * a_f) / 2 -
                    (2 * t3 * j_0) / 3 - (t3 * j_f) / 2;
            nCoeffMat.col(2) =
                    84 * p_0 - 84 * p_f + 45 * t1 * v_0 + 39 * t1 * v_f + 10 * t2 * a_0 - 7 * t2 * a_f + t3 * j_0 +
                    (t3 * j_f) / 2;
            nCoeffMat.col(3) = 35 * p_f - 35 * p_0 - 20 * t1 * v_0 - 15 * t1 * v_f - 5 * t2 * a_0 + (5 * t2 * a_f) / 2 -
                    (2 * t3 * j_0) / 3 - (t3 * j_f) / 6;
            nCoeffMat.col(4) = (t3 * j_0) / 6;
            nCoeffMat.col(5) = (t2 * a_0) / 2;
            nCoeffMat.col(6) = t1 * v_0;
            nCoeffMat.col(7) = p_0;
        } else if (boundCond.cols() == 6) {
            // The BoundaryCond matrix boundCond = [p(0),v(0),a(0),p(T),v(T),a(T)]
            order_ = 5;
            nCoeffMat.resize(3, order_ + 1);
            double t1 = dur;
            double t2 = t1 * t1;

            // Inverse mapping is computed without explicit matrix inverse
            // It maps boundary condition to normalized coefficient matrix
            nCoeffMat.col(0) = 0.5 * (boundCond.col(5) - boundCond.col(2)) * t2 -
                    3.0 * (boundCond.col(1) + boundCond.col(4)) * t1 +
                    6.0 * (boundCond.col(3) - boundCond.col(0));
            nCoeffMat.col(1) = (-boundCond.col(5) + 1.5 * boundCond.col(2)) * t2 +
                    (8.0 * boundCond.col(1) + 7.0 * boundCond.col(4)) * t1 +
                    15.0 * (-boundCond.col(3) + boundCond.col(0));
            nCoeffMat.col(2) = (0.5 * boundCond.col(5) - 1.5 * boundCond.col(2)) * t2 -
                    (6.0 * boundCond.col(1) + 4.0 * boundCond.col(4)) * t1 +
                    10.0 * (boundCond.col(3) - boundCond.col(0));
            nCoeffMat.col(3) = 0.5 * boundCond.col(2) * t2;
            nCoeffMat.col(4) = boundCond.col(1) * t1;
            nCoeffMat.col(5) = boundCond.col(0);
        }
        is_empty_ = false;
    }

    inline bool empty() {
        return is_empty_;
    }

    inline void setCost(double cost) {
        cost_ = cost;
    }

    inline double getCost() {
        return cost_;
    };

    inline void reset() {
        nCoeffMat.setZero();
        duration = 0;
        had_length_ = false;
        had_sampled_ = false;
        sampled_positions_.clear();
    }

    inline void resetDuration(double dur) {
        DynamicMat mat = getCoeffMat();
        duration = dur;
        double t = 1.0;
        cost_ = -1;
        had_sampled_ = false;
        for (int i = order_; i >= 0; i--) {
            nCoeffMat.col(i) = mat.col(i) * t;
            t *= dur;
        }
    }

    inline int getOrder() const {
        return order_;
    }

    inline double getDuration() const {
        return duration;
    }

    // Get the position at time t in this piece
    inline Eigen::Vector3d getPos(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d pos(0.0, 0.0, 0.0);
        double tn = 1.0;
        for (int i = order_; i >= 0; i--) {
            pos += tn * nCoeffMat.col(i);
            tn *= t;
        }
        // The pos is not affected by normalization
        return pos;
    }

    // Get the velocity at time t in this piece
    inline Eigen::Vector3d getVel(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d vel(0.0, 0.0, 0.0);
        double tn = 1.0;
        int n = 1;
        for (int i = order_ - 1; i >= 0; i--) {
            vel += n * tn * nCoeffMat.col(i);
            tn *= t;
            n++;
        }
        // Recover the actual vel
        vel /= duration;
        return vel;
    }

    // Get the acceleration at time t in this piece
    inline Eigen::Vector3d getAcc(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d acc(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        for (int i = order_ - 2; i >= 0; i--) {
            acc += m * n * tn * nCoeffMat.col(i);
            tn *= t;
            m++;
            n++;
        }
        // Recover the actual acc
        acc /= duration * duration;
        return acc;
    }

    // Get the jerk at time t in this piece
    inline Eigen::Vector3d getJerk(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d jerk(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        int k = 3;
        for (int i = order_ - 3; i >= 0; i--) {
            jerk += k * m * n * tn * nCoeffMat.col(i);
            tn *= t;
            k++;
            m++;
            n++;
        }
        // Recover the actual acc
        jerk /= duration * duration * duration;
        return jerk;
    }

    // Get the snap at time t in this piece
    inline Eigen::Vector3d getSnap(double t) const {
        if (t == -1) {
            t = duration;
        }
        // Normalize the time
        t /= duration;
        Eigen::Vector3d snap(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        int k = 3;
        int w = 4;
        for (int i = order_ - 4; i >= 0; i--) {
            snap += w * k * m * n * tn * nCoeffMat.col(i);
            tn *= t;
            w++;
            k++;
            m++;
            n++;
        }
        // Recover the actual acc
        snap /= duration * duration * duration * duration;
        return snap;
    }

    // Get the boundary condition of this piece
    inline DynamicMat getBoundCond() const {
        DynamicMat boundCond;
        if (order_ == 7) {
            boundCond.resize(3, order_ + 1);
            boundCond << getPos(0.0), getVel(0.0), getAcc(0.0), getJerk(0.0),
            getPos(duration), getVel(duration), getAcc(duration), getJerk(duration);
        } else if (order_ == 5) {
            boundCond.resize(3, order_ + 1);
            boundCond << getPos(0.0), getVel(0.0), getAcc(0.0),
            getPos(duration), getVel(duration), getAcc(duration);
        }
        return boundCond;
    }

    // Get the coefficient matrix of the piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline DynamicMat getCoeffMat(bool normalized = false) const {
        DynamicMat posCoeffsMat;
        posCoeffsMat.resize(3, order_ + 1);
        double t = 1;
        for (int i = order_; i >= 0; i--) {
            posCoeffsMat.col(i) = nCoeffMat.col(i) / t;
            t *= normalized ? 1.0 : duration;
        }
        return posCoeffsMat;
    }

    // Get the polynomial coefficients of velocity of this piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline DynamicMat getVelCoeffMat(bool normalized = false) const {
        DynamicMat velCoeffMat;
        velCoeffMat.resize(3, order_);
        int n = 1;
        double t = 1.0;
        t *= normalized ? 1.0 : duration;
        for (int i = order_ - 1; i >= 0; i--) {
            velCoeffMat.col(i) = n * nCoeffMat.col(i) / t;
            n++;
            t *= normalized ? 1.0 : duration;
        }
        return velCoeffMat;
    }

    // Get the polynomial coefficients of acceleration of this piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline DynamicMat getAccCoeffMat(bool normalized = false) const {
        DynamicMat accCoeffMat;
        accCoeffMat.resize(3, order_ - 1);

        int n = 2;
        int m = 1;
        double t = 1.0;
        t *= normalized ? 1.0 : duration * duration;
        for (int i = order_ - 2; i >= 0; i--) {
            accCoeffMat.col(i) = n * m * nCoeffMat.col(i) / t;
            n++;
            m++;
            t *= normalized ? 1.0 : duration;
        }
        return accCoeffMat;
    }

    // Get the max velocity rate of the piece
    inline double getMaxVelRate() const {
        // Compute normalized squared vel norm polynomial coefficient matrix
        Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
        Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                RootFinder::polySqr(nVelCoeffMat.row(1)) +
                RootFinder::polySqr(nVelCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++) {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON) {
            return 0.0;
        } else {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON) {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON) {
                r = 0.5 * (r + 1.0);
            }
            // Find all stationaries
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);

            // Check boundary points and stationaries within duration
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxVelRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
            it != candidates.end();
            it++) {
                if (0.0 <= *it && 1.0 >= *it) {
                    // Recover the actual time then get the vel squared norm
                    tempNormSqr = getVel((*it) * duration).squaredNorm();
                    maxVelRateSqr = maxVelRateSqr < tempNormSqr ? tempNormSqr : maxVelRateSqr;
                }
            }
            return sqrt(maxVelRateSqr);
        }
    }

    // Get the max velocity rate of the piece
    inline double getMinVelRate() const {
        // Compute normalized squared vel norm polynomial coefficient matrix
        Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
        Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                RootFinder::polySqr(nVelCoeffMat.row(1)) +
                RootFinder::polySqr(nVelCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++) {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON) {
            return 0.0;
        } else {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON) {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON) {
                r = 0.5 * (r + 1.0);
            }
            // Find all stationaries
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);

            // Check boundary points and stationaries within duration
            candidates.insert(0.0);
            candidates.insert(1.0);
            double minVelRateSqr = INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
            it != candidates.end();
            it++) {
                if (0.0 <= *it && 1.0 >= *it) {
                    double cur_t = (*it) * duration;
                    if (cur_t < 0.01 || cur_t > duration - 0.01)
                        continue;
                    // Recover the actual time then get the vel squared norm
                    tempNormSqr = getVel(cur_t).squaredNorm();
                    minVelRateSqr = minVelRateSqr > tempNormSqr ? tempNormSqr : minVelRateSqr;
                }
            }
            return sqrt(minVelRateSqr);
        }
    }

    // Get the max acceleration rate of the piece
    inline double getMaxAccRate() const {
        // Compute normalized squared acc norm polynomial coefficient matrix
        Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
        Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                RootFinder::polySqr(nAccCoeffMat.row(1)) +
                RootFinder::polySqr(nAccCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++) {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON) {
            return 0.0;
        } else {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON) {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON) {
                r = 0.5 * (r + 1.0);
            }
            // Find all stationaries
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);
            // Check boundary points and stationaries within duration
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxAccRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
            it != candidates.end();
            it++) {
                if (0.0 <= *it && 1.0 >= *it) {
                    // Recover the actual time then get the acc squared norm
                    tempNormSqr = getAcc((*it) * duration).squaredNorm();
                    maxAccRateSqr = maxAccRateSqr < tempNormSqr ? tempNormSqr : maxAccRateSqr;
                }
            }
            return sqrt(maxAccRateSqr);
        }
    }

    // Check whether velocity rate of the piece is always less than maxVelRate
    inline bool checkMaxVelRate(double maxVelRate) const {
        double sqrMaxVelRate = maxVelRate * maxVelRate;
        if (getVel(0.0).squaredNorm() >= sqrMaxVelRate ||
        getVel(duration).squaredNorm() >= sqrMaxVelRate) {
            return false;
        } else {
            Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
            Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                    RootFinder::polySqr(nVelCoeffMat.row(1)) +
                    RootFinder::polySqr(nVelCoeffMat.row(2));
            // Convert the actual squared maxVelRate to a normalized one
            double t2 = duration * duration;
            coeff.tail<1>()(0) -= sqrMaxVelRate * t2;
            // Directly check the root existence in the normalized interval
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    // Check whether accleration rate of the piece is always less than maxAccRate
    inline bool checkMaxAccRate(double maxAccRate) const {
        double sqrMaxAccRate = maxAccRate * maxAccRate;
        if (getAcc(0.0).squaredNorm() >= sqrMaxAccRate ||
        getAcc(duration).squaredNorm() >= sqrMaxAccRate) {
            return false;
        } else {
            Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
            Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                    RootFinder::polySqr(nAccCoeffMat.row(1)) +
                    RootFinder::polySqr(nAccCoeffMat.row(2));
            // Convert the actual squared maxAccRate to a normalized one
            double t2 = duration * duration;
            double t4 = t2 * t2;
            coeff.tail<1>()(0) -= sqrMaxAccRate * t4;
            // Directly check the root existence in the normalized interval
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    //Scale the Piece(t) to Piece(k*t)
    inline void scaleTime(double k) {
        duration /= k;
        return;
    }

    inline std::vector<Vec3> getTraj(double dt) {
        if (had_sampled_) {
            return sampled_positions_;
        }
        sampled_positions_.clear();
        had_sampled_ = true;
        for (double t = 0.0; t < duration; t += dt) {
            Eigen::Vector3d pos;
            sampled_positions_.push_back(getPos(t));
        }

        return sampled_positions_;
    }

    inline double getLength() {
        if (had_length_) {
            return length_;
        }
        length_ = 0;
        had_length_ = true;
        Vec3 pos = getPos(0), last_pos;
        for (double t = 0.0; t < duration; t += 0.01) {
            last_pos = pos;
            pos = getPos(t);
            length_ += (pos - last_pos).norm();
        }
        return length_;
    }


    /// Return: false means search processes occur error
    enum FW_RETCODE {
        ERROR = -1,
        SUCCESS = 1,
        END_IN_DIS = 2,
        };

    inline int
    getBackwardPosition(const double search_time, ///<[in] input search point, no need on the piece
                        const double dist,    ///<[in
                        double &out_time) {

        /*   /// Typical usage:
         *
             double fw_t;
             int ret_code =  piece.getForwardPosition(search_time, cur_dist, fw_t);
             if(ret_code < 0)
             {
                 print(fg(color::red), " -- [CHK] Fowrdsearch error.\n");
             }
             next_pt = piece.getPos(fw_t);
             search_time = fw_t;
        */
        out_time = -1;
        Vec3 bias_pt = getPos(search_time);
        Eigen::VectorXd coeffsGradT;
        std::set<double> roots;

        if (order_ == 5) {

            coeffsGradT.resize(11);
            DynamicMat coeff_mat = getCoeffMat(), c;
            c.resize(3, 6);

            for (int i = 0; i < 6; i++) {
                c.row(0)[i] = coeff_mat.row(0)[5 - i];
                c.row(1)[i] = coeff_mat.row(1)[5 - i];
                c.row(2)[i] = coeff_mat.row(2)[5 - i];
            }
            c.row(0)[0] -= bias_pt[0];
            c.row(1)[0] -= bias_pt[1];
            c.row(2)[0] -= bias_pt[2];

            coeffsGradT(0) = c(0, 5) * c(0, 5) + c(1, 5) * c(1, 5) + c(2, 5) * c(2, 5);
            coeffsGradT(1) = 2 * c(0, 4) * c(0, 5) + 2 * c(1, 4) * c(1, 5) + 2 * c(2, 4) * c(2, 5);
            coeffsGradT(2) = c(0, 4) * c(0, 4) + c(1, 4) * c(1, 4) + c(2, 4) * c(2, 4) + 2 * c(0, 3) * c(0, 5) +
                    2 * c(1, 3) * c(1, 5) + 2 * c(2, 3) * c(2, 5);
            coeffsGradT(3) =
                    2 * c(0, 2) * c(0, 5) + 2 * c(0, 3) * c(0, 4) + 2 * c(1, 2) * c(1, 5) + 2 * c(1, 3) * c(1, 4) +
                    2 * c(2, 2) * c(2, 5) + 2 * c(2, 3) * c(2, 4);
            coeffsGradT(4) = c(0, 3) * c(0, 3) + c(1, 3) * c(1, 3) + c(2, 3) * c(2, 3) + 2 * c(0, 1) * c(0, 5) +
                    2 * c(0, 2) * c(0, 4) + 2 * c(1, 1) * c(1, 5) + 2 * c(1, 2) * c(1, 4) +
                    2 * c(2, 1) * c(2, 5) + 2 * c(2, 2) * c(2, 4);
            coeffsGradT(5) =
                    2 * c(0, 0) * c(0, 5) + 2 * c(0, 1) * c(0, 4) + 2 * c(0, 2) * c(0, 3) + 2 * c(1, 0) * c(1, 5) +
                    2 * c(1, 1) * c(1, 4) + 2 * c(1, 2) * c(1, 3) + 2 * c(2, 0) * c(2, 5) + 2 * c(2, 1) * c(2, 4) +
                    2 * c(2, 2) * c(2, 3);
            coeffsGradT(6) = c(0, 2) * c(0, 2) + c(1, 2) * c(1, 2) + c(2, 2) * c(2, 2) + 2 * c(0, 0) * c(0, 4) +
                    2 * c(0, 1) * c(0, 3) + 2 * c(1, 0) * c(1, 4) + 2 * c(1, 1) * c(1, 3) +
                    2 * c(2, 0) * c(2, 4) + 2 * c(2, 1) * c(2, 3);
            coeffsGradT(7) =
                    2 * c(0, 0) * c(0, 3) + 2 * c(0, 1) * c(0, 2) + 2 * c(1, 0) * c(1, 3) + 2 * c(1, 1) * c(1, 2) +
                    2 * c(2, 0) * c(2, 3) + 2 * c(2, 1) * c(2, 2);
            coeffsGradT(8) = c(0, 1) * c(0, 1) + c(1, 1) * c(1, 1) + c(2, 1) * c(2, 1) + 2 * c(0, 0) * c(0, 2) +
                    2 * c(1, 0) * c(1, 2) + 2 * c(2, 0) * c(2, 2);
            coeffsGradT(9) = 2 * c(0, 0) * c(0, 1) + 2 * c(1, 0) * c(1, 1) + 2 * c(2, 0) * c(2, 1);
            coeffsGradT(10) = c(0, 0) * c(0, 0) + c(1, 0) * c(1, 0) + c(2, 0) * c(2, 0);
            coeffsGradT(10) -= dist * dist;

            roots = RootFinder::solvePolynomial(coeffsGradT, 0, search_time, 1e-3);

        } else if (order_ == 7) {

            coeffsGradT.resize(15);
            DynamicMat coeff_mat = getCoeffMat(), c;
            c.resize(3, 8);
            Vec3 bias_pt = getPos(search_time);
            for (int i = 0; i < 8; i++) {
                c.row(0)[i] = coeff_mat.row(0)[7 - i];
                c.row(1)[i] = coeff_mat.row(1)[7 - i];
                c.row(2)[i] = coeff_mat.row(2)[7 - i];
            }
            c.row(0)[0] -= bias_pt[0];
            c.row(1)[0] -= bias_pt[1];
            c.row(2)[0] -= bias_pt[2];

            coeffsGradT(0) = c(0, 7) * c(0, 7) + c(1, 7) * c(1, 7) + c(2, 7) * c(2, 7);
            coeffsGradT(1) = 2 * c(0, 6) * c(0, 7) + 2 * c(1, 6) * c(1, 7) + 2 * c(2, 6) * c(2, 7);
            coeffsGradT(2) = c(0, 6) * c(0, 6) + c(1, 6) * c(1, 6) + c(2, 6) * c(2, 6) + 2 * c(0, 5) * c(0, 7) +
                    2 * c(1, 5) * c(1, 7) + 2 * c(2, 5) * c(2, 7);
            coeffsGradT(3) =
                    2 * c(0, 4) * c(0, 7) + 2 * c(0, 5) * c(0, 6) + 2 * c(1, 4) * c(1, 7) + 2 * c(1, 5) * c(1, 6) +
                    2 * c(2, 4) * c(2, 7) + 2 * c(2, 5) * c(2, 6);
            coeffsGradT(4) = c(0, 5) * c(0, 5) + c(1, 5) * c(1, 5) + c(2, 5) * c(2, 5) + 2 * c(0, 3) * c(0, 7) +
                    2 * c(0, 4) * c(0, 6) + 2 * c(1, 3) * c(1, 7) + 2 * c(1, 4) * c(1, 6) +
                    2 * c(2, 3) * c(2, 7) + 2 * c(2, 4) * c(2, 6);
            coeffsGradT(5) =
                    2 * c(0, 2) * c(0, 7) + 2 * c(0, 3) * c(0, 6) + 2 * c(0, 4) * c(0, 5) + 2 * c(1, 2) * c(1, 7) +
                    2 * c(1, 3) * c(1, 6) + 2 * c(1, 4) * c(1, 5) + 2 * c(2, 2) * c(2, 7) + 2 * c(2, 3) * c(2, 6) +
                    2 * c(2, 4) * c(2, 5);
            coeffsGradT(6) = c(0, 4) * c(0, 4) + c(1, 4) * c(1, 4) + c(2, 4) * c(2, 4) + 2 * c(0, 1) * c(0, 7) +
                    2 * c(0, 2) * c(0, 6) + 2 * c(0, 3) * c(0, 5) + 2 * c(1, 1) * c(1, 7) +
                    2 * c(1, 2) * c(1, 6) + 2 * c(1, 3) * c(1, 5) + 2 * c(2, 1) * c(2, 7) +
                    2 * c(2, 2) * c(2, 6) + 2 * c(2, 3) * c(2, 5);
            coeffsGradT(7) =
                    2 * c(0, 0) * c(0, 7) + 2 * c(0, 1) * c(0, 6) + 2 * c(0, 2) * c(0, 5) + 2 * c(0, 3) * c(0, 4) +
                    2 * c(1, 0) * c(1, 7) + 2 * c(1, 1) * c(1, 6) + 2 * c(1, 2) * c(1, 5) + 2 * c(1, 3) * c(1, 4) +
                    2 * c(2, 0) * c(2, 7) + 2 * c(2, 1) * c(2, 6) + 2 * c(2, 2) * c(2, 5) + 2 * c(2, 3) * c(2, 4);
            coeffsGradT(8) = c(0, 3) * c(0, 3) + c(1, 3) * c(1, 3) + c(2, 3) * c(2, 3) + 2 * c(0, 0) * c(0, 6) +
                    2 * c(0, 1) * c(0, 5) + 2 * c(0, 2) * c(0, 4) + 2 * c(1, 0) * c(1, 6) +
                    2 * c(1, 1) * c(1, 5) + 2 * c(1, 2) * c(1, 4) + 2 * c(2, 0) * c(2, 6) +
                    2 * c(2, 1) * c(2, 5) + 2 * c(2, 2) * c(2, 4);
            coeffsGradT(9) =
                    2 * c(0, 0) * c(0, 5) + 2 * c(0, 1) * c(0, 4) + 2 * c(0, 2) * c(0, 3) + 2 * c(1, 0) * c(1, 5) +
                    2 * c(1, 1) * c(1, 4) + 2 * c(1, 2) * c(1, 3) + 2 * c(2, 0) * c(2, 5) + 2 * c(2, 1) * c(2, 4) +
                    2 * c(2, 2) * c(2, 3);
            coeffsGradT(10) = c(0, 2) * c(0, 2) + c(1, 2) * c(1, 2) + c(2, 2) * c(2, 2) + 2 * c(0, 0) * c(0, 4) +
                    2 * c(0, 1) * c(0, 3) + 2 * c(1, 0) * c(1, 4) + 2 * c(1, 1) * c(1, 3) +
                    2 * c(2, 0) * c(2, 4) + 2 * c(2, 1) * c(2, 3);
            coeffsGradT(11) =
                    2 * c(0, 0) * c(0, 3) + 2 * c(0, 1) * c(0, 2) + 2 * c(1, 0) * c(1, 3) + 2 * c(1, 1) * c(1, 2) +
                    2 * c(2, 0) * c(2, 3) + 2 * c(2, 1) * c(2, 2);
            coeffsGradT(12) = c(0, 1) * c(0, 1) + c(1, 1) * c(1, 1) + c(2, 1) * c(2, 1) + 2 * c(0, 0) * c(0, 2) +
                    2 * c(1, 0) * c(1, 2) + 2 * c(2, 0) * c(2, 2);
            coeffsGradT(13) = 2 * c(0, 0) * c(0, 1) + 2 * c(1, 0) * c(1, 1) + 2 * c(2, 0) * c(2, 1);
            coeffsGradT(14) = c(0, 0) * c(0, 0) + c(1, 0) * c(1, 0) + c(2, 0) * c(2, 0);
            coeffsGradT(14) -= dist * dist;
            roots = RootFinder::solvePolynomial(coeffsGradT, 0, search_time, 1e-3);
        } else {
            fmt::print(fg(fmt::color::red), " -- [PIECE] Wrong piece order, force return.\n");
            return ERROR;
        }


        if (roots.size() == 0) {
            fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, " -- [FWS] Cannot get nearest point search.\n");
            return ERROR;
        }

        for (const double &root: roots) {
            if (root >= 0 && root < search_time) {
                out_time = root;
            }
        }
        if (out_time > 0) {
            return SUCCESS;
        }
        fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "All roots larger than duration({}).\n", duration);
        return ERROR;
    }


    inline int
    getNearPosition(const Vec3 search_pt, ///<[in] input search point, no need on the piece
                    const double dist,    ///<[in
                    vector<double> &out_time) {

        /*   /// Typical usage:
         *
             double fw_t;
             int ret_code =  piece.getForwardPosition(search_time, cur_dist, fw_t);
             if(ret_code < 0)
             {
                 print(fg(color::red), " -- [CHK] Fowrdsearch error.\n");
             }
             next_pt = piece.getPos(fw_t);
             search_time = fw_t;
        */

        Vec3 bias_pt = search_pt;
        Eigen::VectorXd coeffsGradT;
        std::set<double> roots;

        if (order_ == 5) {

            coeffsGradT.resize(11);
            DynamicMat coeff_mat = getCoeffMat(), c;
            c.resize(3, 6);

            for (int i = 0; i < 6; i++) {
                c.row(0)[i] = coeff_mat.row(0)[5 - i];
                c.row(1)[i] = coeff_mat.row(1)[5 - i];
                c.row(2)[i] = coeff_mat.row(2)[5 - i];
            }
            c.row(0)[0] -= bias_pt[0];
            c.row(1)[0] -= bias_pt[1];
            c.row(2)[0] -= bias_pt[2];

            coeffsGradT(0) = c(0, 5) * c(0, 5) + c(1, 5) * c(1, 5) + c(2, 5) * c(2, 5);
            coeffsGradT(1) = 2 * c(0, 4) * c(0, 5) + 2 * c(1, 4) * c(1, 5) + 2 * c(2, 4) * c(2, 5);
            coeffsGradT(2) = c(0, 4) * c(0, 4) + c(1, 4) * c(1, 4) + c(2, 4) * c(2, 4) + 2 * c(0, 3) * c(0, 5) +
                    2 * c(1, 3) * c(1, 5) + 2 * c(2, 3) * c(2, 5);
            coeffsGradT(3) =
                    2 * c(0, 2) * c(0, 5) + 2 * c(0, 3) * c(0, 4) + 2 * c(1, 2) * c(1, 5) + 2 * c(1, 3) * c(1, 4) +
                    2 * c(2, 2) * c(2, 5) + 2 * c(2, 3) * c(2, 4);
            coeffsGradT(4) = c(0, 3) * c(0, 3) + c(1, 3) * c(1, 3) + c(2, 3) * c(2, 3) + 2 * c(0, 1) * c(0, 5) +
                    2 * c(0, 2) * c(0, 4) + 2 * c(1, 1) * c(1, 5) + 2 * c(1, 2) * c(1, 4) +
                    2 * c(2, 1) * c(2, 5) + 2 * c(2, 2) * c(2, 4);
            coeffsGradT(5) =
                    2 * c(0, 0) * c(0, 5) + 2 * c(0, 1) * c(0, 4) + 2 * c(0, 2) * c(0, 3) + 2 * c(1, 0) * c(1, 5) +
                    2 * c(1, 1) * c(1, 4) + 2 * c(1, 2) * c(1, 3) + 2 * c(2, 0) * c(2, 5) + 2 * c(2, 1) * c(2, 4) +
                    2 * c(2, 2) * c(2, 3);
            coeffsGradT(6) = c(0, 2) * c(0, 2) + c(1, 2) * c(1, 2) + c(2, 2) * c(2, 2) + 2 * c(0, 0) * c(0, 4) +
                    2 * c(0, 1) * c(0, 3) + 2 * c(1, 0) * c(1, 4) + 2 * c(1, 1) * c(1, 3) +
                    2 * c(2, 0) * c(2, 4) + 2 * c(2, 1) * c(2, 3);
            coeffsGradT(7) =
                    2 * c(0, 0) * c(0, 3) + 2 * c(0, 1) * c(0, 2) + 2 * c(1, 0) * c(1, 3) + 2 * c(1, 1) * c(1, 2) +
                    2 * c(2, 0) * c(2, 3) + 2 * c(2, 1) * c(2, 2);
            coeffsGradT(8) = c(0, 1) * c(0, 1) + c(1, 1) * c(1, 1) + c(2, 1) * c(2, 1) + 2 * c(0, 0) * c(0, 2) +
                    2 * c(1, 0) * c(1, 2) + 2 * c(2, 0) * c(2, 2);
            coeffsGradT(9) = 2 * c(0, 0) * c(0, 1) + 2 * c(1, 0) * c(1, 1) + 2 * c(2, 0) * c(2, 1);
            coeffsGradT(10) = c(0, 0) * c(0, 0) + c(1, 0) * c(1, 0) + c(2, 0) * c(2, 0);
            coeffsGradT(10) -= dist * dist;

            roots = RootFinder::solvePolynomial(coeffsGradT, 0, duration, 1e-3);

        } else if (order_ == 7) {

            coeffsGradT.resize(15);
            DynamicMat coeff_mat = getCoeffMat(), c;
            c.resize(3, 8);
            Vec3 bias_pt = search_pt;
            for (int i = 0; i < 8; i++) {
                c.row(0)[i] = coeff_mat.row(0)[7 - i];
                c.row(1)[i] = coeff_mat.row(1)[7 - i];
                c.row(2)[i] = coeff_mat.row(2)[7 - i];
            }
            c.row(0)[0] -= bias_pt[0];
            c.row(1)[0] -= bias_pt[1];
            c.row(2)[0] -= bias_pt[2];

            coeffsGradT(0) = c(0, 7) * c(0, 7) + c(1, 7) * c(1, 7) + c(2, 7) * c(2, 7);
            coeffsGradT(1) = 2 * c(0, 6) * c(0, 7) + 2 * c(1, 6) * c(1, 7) + 2 * c(2, 6) * c(2, 7);
            coeffsGradT(2) = c(0, 6) * c(0, 6) + c(1, 6) * c(1, 6) + c(2, 6) * c(2, 6) + 2 * c(0, 5) * c(0, 7) +
                    2 * c(1, 5) * c(1, 7) + 2 * c(2, 5) * c(2, 7);
            coeffsGradT(3) =
                    2 * c(0, 4) * c(0, 7) + 2 * c(0, 5) * c(0, 6) + 2 * c(1, 4) * c(1, 7) + 2 * c(1, 5) * c(1, 6) +
                    2 * c(2, 4) * c(2, 7) + 2 * c(2, 5) * c(2, 6);
            coeffsGradT(4) = c(0, 5) * c(0, 5) + c(1, 5) * c(1, 5) + c(2, 5) * c(2, 5) + 2 * c(0, 3) * c(0, 7) +
                    2 * c(0, 4) * c(0, 6) + 2 * c(1, 3) * c(1, 7) + 2 * c(1, 4) * c(1, 6) +
                    2 * c(2, 3) * c(2, 7) + 2 * c(2, 4) * c(2, 6);
            coeffsGradT(5) =
                    2 * c(0, 2) * c(0, 7) + 2 * c(0, 3) * c(0, 6) + 2 * c(0, 4) * c(0, 5) + 2 * c(1, 2) * c(1, 7) +
                    2 * c(1, 3) * c(1, 6) + 2 * c(1, 4) * c(1, 5) + 2 * c(2, 2) * c(2, 7) + 2 * c(2, 3) * c(2, 6) +
                    2 * c(2, 4) * c(2, 5);
            coeffsGradT(6) = c(0, 4) * c(0, 4) + c(1, 4) * c(1, 4) + c(2, 4) * c(2, 4) + 2 * c(0, 1) * c(0, 7) +
                    2 * c(0, 2) * c(0, 6) + 2 * c(0, 3) * c(0, 5) + 2 * c(1, 1) * c(1, 7) +
                    2 * c(1, 2) * c(1, 6) + 2 * c(1, 3) * c(1, 5) + 2 * c(2, 1) * c(2, 7) +
                    2 * c(2, 2) * c(2, 6) + 2 * c(2, 3) * c(2, 5);
            coeffsGradT(7) =
                    2 * c(0, 0) * c(0, 7) + 2 * c(0, 1) * c(0, 6) + 2 * c(0, 2) * c(0, 5) + 2 * c(0, 3) * c(0, 4) +
                    2 * c(1, 0) * c(1, 7) + 2 * c(1, 1) * c(1, 6) + 2 * c(1, 2) * c(1, 5) + 2 * c(1, 3) * c(1, 4) +
                    2 * c(2, 0) * c(2, 7) + 2 * c(2, 1) * c(2, 6) + 2 * c(2, 2) * c(2, 5) + 2 * c(2, 3) * c(2, 4);
            coeffsGradT(8) = c(0, 3) * c(0, 3) + c(1, 3) * c(1, 3) + c(2, 3) * c(2, 3) + 2 * c(0, 0) * c(0, 6) +
                    2 * c(0, 1) * c(0, 5) + 2 * c(0, 2) * c(0, 4) + 2 * c(1, 0) * c(1, 6) +
                    2 * c(1, 1) * c(1, 5) + 2 * c(1, 2) * c(1, 4) + 2 * c(2, 0) * c(2, 6) +
                    2 * c(2, 1) * c(2, 5) + 2 * c(2, 2) * c(2, 4);
            coeffsGradT(9) =
                    2 * c(0, 0) * c(0, 5) + 2 * c(0, 1) * c(0, 4) + 2 * c(0, 2) * c(0, 3) + 2 * c(1, 0) * c(1, 5) +
                    2 * c(1, 1) * c(1, 4) + 2 * c(1, 2) * c(1, 3) + 2 * c(2, 0) * c(2, 5) + 2 * c(2, 1) * c(2, 4) +
                    2 * c(2, 2) * c(2, 3);
            coeffsGradT(10) = c(0, 2) * c(0, 2) + c(1, 2) * c(1, 2) + c(2, 2) * c(2, 2) + 2 * c(0, 0) * c(0, 4) +
                    2 * c(0, 1) * c(0, 3) + 2 * c(1, 0) * c(1, 4) + 2 * c(1, 1) * c(1, 3) +
                    2 * c(2, 0) * c(2, 4) + 2 * c(2, 1) * c(2, 3);
            coeffsGradT(11) =
                    2 * c(0, 0) * c(0, 3) + 2 * c(0, 1) * c(0, 2) + 2 * c(1, 0) * c(1, 3) + 2 * c(1, 1) * c(1, 2) +
                    2 * c(2, 0) * c(2, 3) + 2 * c(2, 1) * c(2, 2);
            coeffsGradT(12) = c(0, 1) * c(0, 1) + c(1, 1) * c(1, 1) + c(2, 1) * c(2, 1) + 2 * c(0, 0) * c(0, 2) +
                    2 * c(1, 0) * c(1, 2) + 2 * c(2, 0) * c(2, 2);
            coeffsGradT(13) = 2 * c(0, 0) * c(0, 1) + 2 * c(1, 0) * c(1, 1) + 2 * c(2, 0) * c(2, 1);
            coeffsGradT(14) = c(0, 0) * c(0, 0) + c(1, 0) * c(1, 0) + c(2, 0) * c(2, 0);
            coeffsGradT(14) -= dist * dist;
            roots = RootFinder::solvePolynomial(coeffsGradT, 0, duration, 1e-3);
        } else {
            fmt::print(fg(fmt::color::red), " -- [PIECE] Wrong piece order, force return.\n");
            return ERROR;
        }


        if (roots.size() == 0) {
            out_time.clear();
            fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, " -- [FWS] Cannot get nearest point search.\n");
            return ERROR;
        }

        for (const double &root: roots) {
            if (root >= 0 && root < duration) {
                out_time.push_back(root);
            }
        }
        if (out_time.size() > 0) {
            return SUCCESS;
        }
        fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "All roots larger than duration({}).\n", duration);
        return ERROR;
    }


    /// Do forward search in [search_time, duration] all t out of this range is blocked.
    inline int
    getForwardPosition(const double search_time, ///<[in] input search time on the piece, should be [0, duration]'
                       const double dist,    ///<[in
                       double &out_time) {

        /*   /// Typical usage:
         *
             double fw_t;
             int ret_code =  piece.getForwardPosition(search_time, cur_dist, fw_t);
             if(ret_code < 0)
             {
                 print(fg(color::red), " -- [CHK] Fowrdsearch error.\n");
             }
             next_pt = piece.getPos(fw_t);
             search_time = fw_t;
        */

        if (search_time > duration) {
            fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, " -- [FWS] Forward search time out of traj.\n");
            return ERROR;
        }

        Vec3 bias_pt = getPos(search_time);
        Eigen::VectorXd coeffsGradT;
        std::set<double> roots;

        if (order_ == 5) {

            coeffsGradT.resize(11);
            DynamicMat coeff_mat = getCoeffMat(), c;
            c.resize(3, 6);

            for (int i = 0; i < 6; i++) {
                c.row(0)[i] = coeff_mat.row(0)[5 - i];
                c.row(1)[i] = coeff_mat.row(1)[5 - i];
                c.row(2)[i] = coeff_mat.row(2)[5 - i];
            }
            c.row(0)[0] -= bias_pt[0];
            c.row(1)[0] -= bias_pt[1];
            c.row(2)[0] -= bias_pt[2];

            coeffsGradT(0) = c(0, 5) * c(0, 5) + c(1, 5) * c(1, 5) + c(2, 5) * c(2, 5);
            coeffsGradT(1) = 2 * c(0, 4) * c(0, 5) + 2 * c(1, 4) * c(1, 5) + 2 * c(2, 4) * c(2, 5);
            coeffsGradT(2) = c(0, 4) * c(0, 4) + c(1, 4) * c(1, 4) + c(2, 4) * c(2, 4) + 2 * c(0, 3) * c(0, 5) +
                    2 * c(1, 3) * c(1, 5) + 2 * c(2, 3) * c(2, 5);
            coeffsGradT(3) =
                    2 * c(0, 2) * c(0, 5) + 2 * c(0, 3) * c(0, 4) + 2 * c(1, 2) * c(1, 5) + 2 * c(1, 3) * c(1, 4) +
                    2 * c(2, 2) * c(2, 5) + 2 * c(2, 3) * c(2, 4);
            coeffsGradT(4) = c(0, 3) * c(0, 3) + c(1, 3) * c(1, 3) + c(2, 3) * c(2, 3) + 2 * c(0, 1) * c(0, 5) +
                    2 * c(0, 2) * c(0, 4) + 2 * c(1, 1) * c(1, 5) + 2 * c(1, 2) * c(1, 4) +
                    2 * c(2, 1) * c(2, 5) + 2 * c(2, 2) * c(2, 4);
            coeffsGradT(5) =
                    2 * c(0, 0) * c(0, 5) + 2 * c(0, 1) * c(0, 4) + 2 * c(0, 2) * c(0, 3) + 2 * c(1, 0) * c(1, 5) +
                    2 * c(1, 1) * c(1, 4) + 2 * c(1, 2) * c(1, 3) + 2 * c(2, 0) * c(2, 5) + 2 * c(2, 1) * c(2, 4) +
                    2 * c(2, 2) * c(2, 3);
            coeffsGradT(6) = c(0, 2) * c(0, 2) + c(1, 2) * c(1, 2) + c(2, 2) * c(2, 2) + 2 * c(0, 0) * c(0, 4) +
                    2 * c(0, 1) * c(0, 3) + 2 * c(1, 0) * c(1, 4) + 2 * c(1, 1) * c(1, 3) +
                    2 * c(2, 0) * c(2, 4) + 2 * c(2, 1) * c(2, 3);
            coeffsGradT(7) =
                    2 * c(0, 0) * c(0, 3) + 2 * c(0, 1) * c(0, 2) + 2 * c(1, 0) * c(1, 3) + 2 * c(1, 1) * c(1, 2) +
                    2 * c(2, 0) * c(2, 3) + 2 * c(2, 1) * c(2, 2);
            coeffsGradT(8) = c(0, 1) * c(0, 1) + c(1, 1) * c(1, 1) + c(2, 1) * c(2, 1) + 2 * c(0, 0) * c(0, 2) +
                    2 * c(1, 0) * c(1, 2) + 2 * c(2, 0) * c(2, 2);
            coeffsGradT(9) = 2 * c(0, 0) * c(0, 1) + 2 * c(1, 0) * c(1, 1) + 2 * c(2, 0) * c(2, 1);
            coeffsGradT(10) = c(0, 0) * c(0, 0) + c(1, 0) * c(1, 0) + c(2, 0) * c(2, 0);
            coeffsGradT(10) -= dist * dist;

            roots = RootFinder::solvePolynomial(coeffsGradT, search_time, duration, 1e-3);

        } else if (order_ == 7) {

            coeffsGradT.resize(15);
            DynamicMat coeff_mat = getCoeffMat(), c;
            c.resize(3, 8);
            Vec3 bias_pt = getPos(search_time);
            for (int i = 0; i < 8; i++) {
                c.row(0)[i] = coeff_mat.row(0)[7 - i];
                c.row(1)[i] = coeff_mat.row(1)[7 - i];
                c.row(2)[i] = coeff_mat.row(2)[7 - i];
            }
            c.row(0)[0] -= bias_pt[0];
            c.row(1)[0] -= bias_pt[1];
            c.row(2)[0] -= bias_pt[2];

            coeffsGradT(0) = c(0, 7) * c(0, 7) + c(1, 7) * c(1, 7) + c(2, 7) * c(2, 7);
            coeffsGradT(1) = 2 * c(0, 6) * c(0, 7) + 2 * c(1, 6) * c(1, 7) + 2 * c(2, 6) * c(2, 7);
            coeffsGradT(2) = c(0, 6) * c(0, 6) + c(1, 6) * c(1, 6) + c(2, 6) * c(2, 6) + 2 * c(0, 5) * c(0, 7) +
                    2 * c(1, 5) * c(1, 7) + 2 * c(2, 5) * c(2, 7);
            coeffsGradT(3) =
                    2 * c(0, 4) * c(0, 7) + 2 * c(0, 5) * c(0, 6) + 2 * c(1, 4) * c(1, 7) + 2 * c(1, 5) * c(1, 6) +
                    2 * c(2, 4) * c(2, 7) + 2 * c(2, 5) * c(2, 6);
            coeffsGradT(4) = c(0, 5) * c(0, 5) + c(1, 5) * c(1, 5) + c(2, 5) * c(2, 5) + 2 * c(0, 3) * c(0, 7) +
                    2 * c(0, 4) * c(0, 6) + 2 * c(1, 3) * c(1, 7) + 2 * c(1, 4) * c(1, 6) +
                    2 * c(2, 3) * c(2, 7) + 2 * c(2, 4) * c(2, 6);
            coeffsGradT(5) =
                    2 * c(0, 2) * c(0, 7) + 2 * c(0, 3) * c(0, 6) + 2 * c(0, 4) * c(0, 5) + 2 * c(1, 2) * c(1, 7) +
                    2 * c(1, 3) * c(1, 6) + 2 * c(1, 4) * c(1, 5) + 2 * c(2, 2) * c(2, 7) + 2 * c(2, 3) * c(2, 6) +
                    2 * c(2, 4) * c(2, 5);
            coeffsGradT(6) = c(0, 4) * c(0, 4) + c(1, 4) * c(1, 4) + c(2, 4) * c(2, 4) + 2 * c(0, 1) * c(0, 7) +
                    2 * c(0, 2) * c(0, 6) + 2 * c(0, 3) * c(0, 5) + 2 * c(1, 1) * c(1, 7) +
                    2 * c(1, 2) * c(1, 6) + 2 * c(1, 3) * c(1, 5) + 2 * c(2, 1) * c(2, 7) +
                    2 * c(2, 2) * c(2, 6) + 2 * c(2, 3) * c(2, 5);
            coeffsGradT(7) =
                    2 * c(0, 0) * c(0, 7) + 2 * c(0, 1) * c(0, 6) + 2 * c(0, 2) * c(0, 5) + 2 * c(0, 3) * c(0, 4) +
                    2 * c(1, 0) * c(1, 7) + 2 * c(1, 1) * c(1, 6) + 2 * c(1, 2) * c(1, 5) + 2 * c(1, 3) * c(1, 4) +
                    2 * c(2, 0) * c(2, 7) + 2 * c(2, 1) * c(2, 6) + 2 * c(2, 2) * c(2, 5) + 2 * c(2, 3) * c(2, 4);
            coeffsGradT(8) = c(0, 3) * c(0, 3) + c(1, 3) * c(1, 3) + c(2, 3) * c(2, 3) + 2 * c(0, 0) * c(0, 6) +
                    2 * c(0, 1) * c(0, 5) + 2 * c(0, 2) * c(0, 4) + 2 * c(1, 0) * c(1, 6) +
                    2 * c(1, 1) * c(1, 5) + 2 * c(1, 2) * c(1, 4) + 2 * c(2, 0) * c(2, 6) +
                    2 * c(2, 1) * c(2, 5) + 2 * c(2, 2) * c(2, 4);
            coeffsGradT(9) =
                    2 * c(0, 0) * c(0, 5) + 2 * c(0, 1) * c(0, 4) + 2 * c(0, 2) * c(0, 3) + 2 * c(1, 0) * c(1, 5) +
                    2 * c(1, 1) * c(1, 4) + 2 * c(1, 2) * c(1, 3) + 2 * c(2, 0) * c(2, 5) + 2 * c(2, 1) * c(2, 4) +
                    2 * c(2, 2) * c(2, 3);
            coeffsGradT(10) = c(0, 2) * c(0, 2) + c(1, 2) * c(1, 2) + c(2, 2) * c(2, 2) + 2 * c(0, 0) * c(0, 4) +
                    2 * c(0, 1) * c(0, 3) + 2 * c(1, 0) * c(1, 4) + 2 * c(1, 1) * c(1, 3) +
                    2 * c(2, 0) * c(2, 4) + 2 * c(2, 1) * c(2, 3);
            coeffsGradT(11) =
                    2 * c(0, 0) * c(0, 3) + 2 * c(0, 1) * c(0, 2) + 2 * c(1, 0) * c(1, 3) + 2 * c(1, 1) * c(1, 2) +
                    2 * c(2, 0) * c(2, 3) + 2 * c(2, 1) * c(2, 2);
            coeffsGradT(12) = c(0, 1) * c(0, 1) + c(1, 1) * c(1, 1) + c(2, 1) * c(2, 1) + 2 * c(0, 0) * c(0, 2) +
                    2 * c(1, 0) * c(1, 2) + 2 * c(2, 0) * c(2, 2);
            coeffsGradT(13) = 2 * c(0, 0) * c(0, 1) + 2 * c(1, 0) * c(1, 1) + 2 * c(2, 0) * c(2, 1);
            coeffsGradT(14) = c(0, 0) * c(0, 0) + c(1, 0) * c(1, 0) + c(2, 0) * c(2, 0);
            coeffsGradT(14) -= dist * dist;
            roots = RootFinder::solvePolynomial(coeffsGradT, search_time, duration, 1e-3);
        } else {
            fmt::print(fg(fmt::color::red), " -- [PIECE] Wrong piece order, force return.\n");
            return ERROR;
        }


        if (roots.size() == 0) {
            double end_dis = (bias_pt - getPos(-1)).norm();
            if (end_dis < dist + 0.01) {
                out_time = duration;
                return END_IN_DIS;
            }
            fmt::print(" -- [FWS] Cur t = {}, total_t = {} , forward_dis = {}, end_dis = {}\n", search_time, duration,
                       dist,
                       end_dis);
            fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, " -- [FWS] Cannot get forward search.\n");
            return ERROR;
        }

        for (const double &root: roots) {
            if (root > duration) {
                continue;
            }
            if (root - search_time >= 0) {
                out_time = root;
                return SUCCESS;
            }
        }
        fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "All roots larger than duration({}).\n", duration);
        return ERROR;
    }

};

#endif //POLY_DATA_STRUCTURE_PIECE_HPP
