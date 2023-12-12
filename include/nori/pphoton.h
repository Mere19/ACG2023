#include <nori/kdtree.h>
#include <nori/color.h>

/// progressive photon class

NORI_NAMESPACE_BEGIN

struct PPhotonData {
    Vector3f n;     // normal at position x
    Vector3f dir;   // ray direction
    float R;        // current photon radius
    int N;          // accumulated photon count
    Color3f power;  // accumulated power

    PPhotonData() { }

    PPhotonData(const Vector3f &n, const Vector3f &dir, const Color3f &power);
};

struct PPhoton : public GenericKDTreeNode<Point3f, PPhotonData> {
public:
    /// dummy constructor
    PPhoton() { }

    /// construct from hitpoint
    PPhoton(const Point3f &pos, const Vector3f &n, const Vector3f &dir, const Color3f &power)
        : GenericKDTreeNode(pos, PPhotonData(n, dir, power)) { }
};

NORI_NAMESPACE_END