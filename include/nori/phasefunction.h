#if !defined(__NORI_PF_H)
#define __NORI_PF_H
#include <nori/object.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN
struct PhaseFunctionQueryRecord {
    Vector3f wi;
    Vector3f wo;
    float pdf;
    EMeasure measure;
    PhaseFunctionQueryRecord() {}
    PhaseFunctionQueryRecord(const Vector3f& wi)
        : wi(wi), measure(EUnknownMeasure) { }
    PhaseFunctionQueryRecord(const Vector3f& wi,
        const Vector3f& wo, EMeasure measure)
        : wi(wi), wo(wo), measure(measure) { }
};

class PhaseFunction : public NoriObject {
public:
    virtual float sample(PhaseFunctionQueryRecord& pRec, const Point2f& sample) const = 0;
    virtual float pdf(PhaseFunctionQueryRecord& pRec) const = 0;
    virtual float eval(PhaseFunctionQueryRecord& pRec) const = 0;
    EClassType getClassType() const { return EPhaseFunction; }

protected:
};

NORI_NAMESPACE_END
#endif