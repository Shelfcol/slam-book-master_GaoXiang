/**
 * @file GPSPose2Factor.h
 * @brief 2D 'GPS' like factor for Pose2
 * @date Nov 8, 2016
 * @author Jing Dong
 */

/**
 * A simple 2D 'GPS' like factor
 * The factor contains a X-Y position measurement (mx, my) for a Pose2, but no rotation information
 * The error vector will be [x-mx, y-my]'
 */

#pragma once

#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/base/Matrix.h>
#include <gtsam/base/Vector.h>
#include <gtsam/geometry/Pose2.h>


// you can custom namespace (if needed for your project)
namespace gtsamexamples {


class GPSPose2Factor: public gtsam::NoiseModelFactor1<gtsam::Pose2> {//Factor1表示只会连接一个factor，因为gtsam中里程计的值和landmark的观测值才是中间的小黑点，而机器人状态和landmark的坐标是factor。这里GPS的测量值只是小黑点

private:
  // measurement information
  double mx_, my_;

public:

  /**
   * Constructor
   * @param poseKey    associated pose varible key
   * @param model      noise model for GPS snesor, in X-Y
   * @param m          Point2 measurement
   */
  GPSPose2Factor(gtsam::Key poseKey, const gtsam::Point2 m, gtsam::SharedNoiseModel model) :
      gtsam::NoiseModelFactor1<gtsam::Pose2>(model, poseKey), mx_(m.x()), my_(m.y()) {}

  // error function
  // @param p    the pose in Pose2
  // @param H    the optional Jacobian matrix, which use boost optional and has default null pointer
  gtsam::Vector evaluateError(const gtsam::Pose2& p, boost::optional<gtsam::Matrix&> H = boost::none) const {
  
    // note that use boost optional like a pointer
    // only calculate jacobian matrix when non-null pointer exists
    if (H) *H = (gtsam::Matrix23() << 1.0, 0.0, 0.0, //Jacobian矩阵！！！！！！！！！！！
                                      0.0, 1.0, 0.0).finished();
    
    // return error vector
    return (gtsam::Vector2() << p.x() - mx_, p.y() - my_).finished();
  }

};

} // namespace gtsamexamples

