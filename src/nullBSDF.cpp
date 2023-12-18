#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

class NullBSDF : public BSDF {
public:
    NullBSDF(const PropertyList&) { }

    virtual Color3f eval(const BSDFQueryRecord&) const override {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    virtual float pdf(const BSDFQueryRecord&) const override {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    virtual Color3f sample(BSDFQueryRecord& bRec, const Point2f&) const override {

        // Reflection in local coordinates
        bRec.wo = Vector3f(
            -bRec.wi.x(),
            -bRec.wi.y(),
            -bRec.wi.z()
        );
        bRec.measure = EDiscrete;

        return Color3f(1.0f);
    }

    virtual std::string toString() const override {
        return "NULLBSDF[]";
    }
};

NORI_REGISTER_CLASS(NullBSDF, "nullBSDF");
NORI_NAMESPACE_END
