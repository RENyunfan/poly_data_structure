//
// Created by yunfan on 2021/11/6.
//

#ifndef POLY_DATA_STRUCTURE_TRAJECTORY_HPP
#define POLY_DATA_STRUCTURE_TRAJECTORY_HPP

#include "common_include.hpp"
#include "Piece.hpp"
#include "BandedSystem.hpp"

// A whole trajectory which contains multiple pieces
class Trajectory {
private:
    typedef std::vector<Piece> Pieces;
    Pieces pieces;
    bool had_sampled_ = false;
    std::vector<Vec3> sampled_positions_;
    std::vector<Ball> sfcs;
    bool had_length_ = false;
    double length_;
    double cost_;

public:
    Trajectory() = default;

    // Constructor from durations and coefficient matrices
    Trajectory(const std::vector<double> &durs,
               const std::vector<DynamicMat> &coeffMats) {
        int N = std::min(durs.size(), coeffMats.size());
        for (int i = 0; i < N; i++) {
            pieces.emplace_back(durs[i], coeffMats[i]);
        }
    }

    inline bool empty(){
        return getTotalDuration() > 0.01?false:true;
    }

    inline void setCost(double cost) {
        cost_ = cost;
    }

    inline double getCost() {
        return cost_;
    }

    inline size_t getPieceNum() const {
        return pieces.size();
    }

    inline void reset() {
        pieces.clear();
        had_length_ = false;
        had_sampled_ = false;
        sampled_positions_.clear();
    }

    // Get durations vector of all pieces
    inline std::vector<double> getDurations() const {
        std::vector<double> durations;
        durations.reserve(getPieceNum());
        for (int i = 0; i < getPieceNum(); i++) {
            durations.push_back(pieces[i].getDuration());
        }
        return durations;
    }

    // Get total duration of the trajectory
    inline double getTotalDuration() const {
        double totalDuration = 0.0;
        for (int i = 0; i < getPieceNum(); i++) {
            totalDuration += pieces[i].getDuration();
        }
        return totalDuration;
    }

    inline double getTotalLength() {
        if (had_length_) {
            return length_;
        }
        length_ = 0;
        had_length_ = true;
        int num_seg = getPieceNum();
        for (size_t t = 0; t < num_seg; t++) {
            length_ += pieces[t].getLength();
        }
        return length_;
    }

    // Reload the operator[] to reach the i-th piece
    inline const Piece &operator[](int i) const {
        return pieces[i];
    }

    inline Piece &operator[](int i) {
        return pieces[i];
    }

    inline void clear(void) {
        had_length_ = false;
        pieces.clear();
    }

    inline Pieces::const_iterator begin() const {
        return pieces.begin();
    }

    inline Pieces::const_iterator end() const {
        return pieces.end();
    }

    // Put another piece at the tail of this trajectory
    inline void emplace_back(const Piece &piece) {
        had_length_ = false;
        pieces.emplace_back(piece);
        return;
    }

    // Two corresponding constructors of Piece both are supported here
    template<typename ArgTypeL, typename ArgTypeR>
    inline void emplace_back(const ArgTypeL &argL, const ArgTypeR &argR) {
        had_length_ = false;
        pieces.emplace_back(argL, argR);
        return;
    }

    // Append another Trajectory at the tail of this trajectory
    inline void append(const Trajectory &traj) {
        had_length_ = false;
        pieces.insert(pieces.end(), traj.begin(), traj.end());
        return;
    }

    inline double getSegmentDuration(int start_id, int end_id){
        double seg_dur = 0.0;

        for(int i = start_id; i <=end_id;i++){
            seg_dur += pieces[i].getDuration();
        }
        return seg_dur;
    }

    // Find the piece at which the time t is located
    // The index is returned and the offset in t is removed
    inline int locatePieceIdx(double &t) const {
        int idx;
        double dur;
        for (idx = 0;
        idx < getPieceNum() &&
        t > (dur = pieces[idx].getDuration());
        idx++) {
            t -= dur;
        }
        if (idx == getPieceNum()) {
            idx--;
            t += pieces[idx].getDuration();
        }
        return idx;
    }

    // Get the position at time t of the trajectory
    inline Eigen::Vector3d getPos(double t) const {
        if (t == -1) {
            return pieces[getPieceNum() - 1].getPos(-1);
        }
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getPos(t);
    }

    // Get the velocity at time t of the trajectory
    inline Eigen::Vector3d getVel(double t) const {
        if (t == -1) {
            return pieces[getPieceNum() - 1].getVel(-1);
        }
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getVel(t);
    }

    // Get the acceleration at time t of the trajectory
    inline Eigen::Vector3d getAcc(double t) const {
        if (t == -1) {
            return pieces[getPieceNum() - 1].getAcc(-1);
        }
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getAcc(t);
    }

    // Get the acceleration at time t of the trajectory
    inline Eigen::Vector3d getJerk(double t) const {
        if (t == -1) {
            return pieces[getPieceNum() - 1].getJerk(-1);
        }
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getJerk(t);
    }

    // Get the acceleration at time t of the trajectory
    inline Eigen::Vector3d getSnap(double t) const {
        if (t == -1) {
            return pieces[getPieceNum() - 1].getSnap(-1);
        }
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getSnap(t);
    }

    inline Eigen::MatrixXd getState(double t) const {
        StatePVAJ state;
        state << getPos(t), getVel(t), getAcc(t), getJerk(t);
        return state;
    }

    // Get the position at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncPos(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getPos(0.0);
        } else {
            return pieces[juncIdx - 1].getPos(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the velocity at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncVel(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getVel(0.0);
        } else {
            return pieces[juncIdx - 1].getVel(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the acceleration at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncAcc(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getAcc(0.0);
        } else {
            return pieces[juncIdx - 1].getAcc(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the acceleration at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncJerk(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getJerk(0.0);
        } else {
            return pieces[juncIdx - 1].getJerk(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the max velocity rate of the trajectory
    inline double getMaxVelRate() const {
        double maxVelRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++) {
            tempNorm = pieces[i].getMaxVelRate();
            maxVelRate = maxVelRate < tempNorm ? tempNorm : maxVelRate;
        }
        return maxVelRate;
    }

    // Get the min velocity rate of the trajectory
    inline double getMinVelRate() const {
        double minVelRate = INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++) {
            tempNorm = pieces[i].getMinVelRate();
            minVelRate = minVelRate > tempNorm ? tempNorm : minVelRate;
        }
        return minVelRate;
    }

    // Get the max acceleration rate of the trajectory
    inline double getMaxAccRate() const {
        double maxAccRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++) {
            tempNorm = pieces[i].getMaxAccRate();
            maxAccRate = maxAccRate < tempNorm ? tempNorm : maxAccRate;
        }
        return maxAccRate;
    }

    // Check whether the velocity rate of this trajectory exceeds the threshold
    inline bool checkMaxVelRate(double maxVelRate) const {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++) {
            feasible = feasible && pieces[i].checkMaxVelRate(maxVelRate);
        }
        return feasible;
    }

    // Check whether the acceleration rate of this trajectory exceeds the threshold
    inline bool checkMaxAccRate(double maxAccRate) const {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++) {
            feasible = feasible && pieces[i].checkMaxAccRate(maxAccRate);
        }
        return feasible;
    }

    // Scale the Trajectory(t) to Trajectory(k*t)
    inline void scaleTime(double k) {
        for (int i = 0; i < getPieceNum(); i++) {
            pieces[i].scaleTime(k);
        }
    }

    inline void reserve(const int &n) {
        pieces.reserve(n);
        return;
    }


    /* For optimization Usage */

    inline std::vector<Eigen::Vector3d> getUniformWaypoints(int num_per_piece = 2) {
        int num_seg = getPieceNum();
        int total_num = num_per_piece * num_seg + num_seg - 1;
        double dur = getTotalDuration();
        double seg_t = dur / total_num;
        double eval_t = seg_t;
        vector<Vec3> waypts;
        while (eval_t < dur) {
            waypts.push_back(getPos(eval_t));
            eval_t += seg_t;
            if (dur - eval_t < seg_t) {
                break;
            }
        }
        return waypts;
    }

    inline std::vector<Eigen::Vector3d> getWaypoints() {

        std::vector<Eigen::Vector3d> waypts;
        int num_seg = getPieceNum();
        for (size_t i = 1; i < num_seg; i++) {
            waypts.push_back(pieces[i].getPos(0));
        }

        /* double cur_dur, total_dur = getTotalDuration();
         for (size_t i = 0; i < num_seg; i++) {
             waypts.push_back(pieces[i].getPos(0));
             if (!pure_waypts) {
                 cur_dur = pieces[i].getDuration();
                 if (cur_dur > 0.4 * total_dur) {
                     waypts.push_back(pieces[i].getPos(cur_dur / 3));
                     waypts.push_back(pieces[i].getPos(cur_dur / 3 * 2));
                 }
             }
         }
         waypts.push_back(getPos(total_dur));*/
        return waypts;
    }



    /* For UAV Usage */

    inline void getRotation(const double &t,
                            const double &yaw,
                            const double &gAcc,
                            Eigen::Matrix3d &rotM) const {
        rotM.col(2) = getAcc(t);
        rotM(2, 2) += gAcc;
        rotM.col(2).normalize();
        rotM.col(1) = rotM.col(2).cross(Eigen::Vector3d(cos(yaw), sin(yaw), 0.0));
        rotM.col(1).normalize();
        rotM.col(0) = rotM.col(1).cross(rotM.col(2));
        return;
    }

    inline void getOmegaAndAtAndRpy(const double &t,
                                    const double &yaw,
                                    const double &yaw_dot,
                                    const double &gAcc,
                                    Vec3 &omega,
                                    Vec3 &rpy,
                                    double &aT) const {
        Vec3 cur_jerk = getJerk(t);
        Vec3 gravity(0, 0, gAcc);
        Vec3 cur_acc = getAcc(t);
        aT = (gravity + cur_acc).norm();


        Vec3 xB, yB, zB;
        Vec3 xC(cos(yaw), sin(yaw), 0);

        zB = (gravity + cur_acc).normalized();
        yB = ((zB).cross(xC)).normalized();
        xB = yB.cross(zB);
        Eigen::Matrix<double, 3, 3> rotation_matrix;
        rotation_matrix << xB, yB, zB;
        rpy = rotation_matrix.eulerAngles(2, 1, 0);
        Vec3 hw = (cur_jerk - (zB.dot(cur_jerk) * zB)) / aT;
        omega(0) = hw.dot(yB);
        omega(1) = hw.dot(xB);
        omega(2) = yaw_dot * (zB.dot(Vec3(0, 0, 1)));

    }

    inline Eigen::Vector3d getTiltRate(const double &t,
                                       const double &gAcc) const {
        Eigen::Matrix3d rotM;
        rotM.col(2) = getAcc(t);
        rotM(2, 2) += gAcc;
        double thrAcc = rotM.col(2).norm();
        rotM.col(2) /= thrAcc;
        rotM.col(1) = rotM.col(2).cross(Eigen::Vector3d::UnitX());
        rotM.col(1).normalize();
        rotM.col(0) = rotM.col(1).cross(rotM.col(2));
        Eigen::Vector3d bdr = rotM.transpose() * getJerk(t) / thrAcc;
        return Eigen::Vector3d(-bdr(1), bdr(0), 0.0);
    }
};


#endif //POLY_DATA_STRUCTURE_TRAJECTORY_HPP
