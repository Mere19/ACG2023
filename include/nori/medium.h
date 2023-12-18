#if !defined(__NORI_MEDIUM_H)
#define __NORI_MEDIUM_H
#include <nori/phasefunction.h>
#include <nori/object.h>
class Shape;
class Volume;

NORI_NAMESPACE_BEGIN
    struct MediumQueryRecord {
    // reference point
    Point3f ref;
    /// sampled point
    Point3f p;
    // incoming Direction
    Vector3f wi;
    // Sampled direction
    Vector3f wo;
    // Max free path length
    float tMax;
    // return value based on the queried point
    Color3f ret;

    //illumination on the sampled point
    Color3f radiance;

    float radiance_pdf;

    EMeasure measure;

    Ray3f shadowRay;

    bool debug;

    /// Empty constructor
    MediumQueryRecord() {}
    /// Sample phase function and distance (free path)
    MediumQueryRecord(const Point3f& _ref, const Vector3f& _wi, const float& _tMax) : ref(_ref), wi(_wi), tMax(_tMax) {
        debug = false;
    }
    /// Evaluate albedo
    MediumQueryRecord(const Point3f& _ref) : ref(_ref) {}
    /// Evaluate transmittance
    MediumQueryRecord(const Point3f& _ref, const Point3f& _p) : ref(_ref), p(_p) {
        wo = (p - ref).normalized();
    }
    /// Query probability density of sampling
    MediumQueryRecord(const Point3f& _ref, const Point3f& _p, const Vector3f& _wi, const EMeasure& _measure) : ref(_ref), p(_p), wi(_wi), measure(_measure) {
        wo = (p - ref).normalized();
    }
};

/**
 * \brief Superclass of all mediums
 */
class Medium : public NoriObject {
public:
    virtual  bool sample_freepath(MediumQueryRecord& mRec, Sampler* sampler) const = 0;

    virtual Color3f Tr(const MediumQueryRecord& mRec, Sampler* sampler) const { return Color3f(0); }

    // Return the phase function of this medium
    inline const PhaseFunction* getPhaseFunction() const { return m_phase; }

    //property check
    virtual bool isHomogeneous() const { return false; }
    virtual bool isHeterogeneous() const { return false; }
    virtual bool isEmissive() const { return false; }

    //eval
    virtual Color3f eval_radiance(Point3f& p) const { return Color3f(0); }

    //sample
    virtual Color3f sample_radiance(MediumQueryRecord& mRec, Sampler* sampler) const { return Color3f(0); }

    virtual void setShape(const Shape* shape) { m_shape = shape; }
    virtual const Shape* getShape() const { return m_shape; }

    virtual void addChild(NoriObject* obj) {}

    virtual std::string toString() const = 0;

    //generate the volume grid for perlin noise
    virtual void volGrid() {};

    EClassType getClassType() const override { return EMedium; }



protected:
    PhaseFunction* m_phase = nullptr;
    const Shape* m_shape = nullptr;
    // Color3f m_sigmaA;
    // Color3f m_sigmaS;
    // Color3f m_sigmaT;
};

NORI_NAMESPACE_END
#endif /* __NORI_MEDIUM_H */