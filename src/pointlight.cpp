#include<nori/emitter.h>
#include<nori/sampler.h>
#include<nori/ray.h>

NORI_NAMESPACE_BEGIN

class PointLight : public Emitter{
private:
    /* data */
    Color3f power;
    Point3f position;
public:
    PointLight(const PropertyList & propList) {
        power = propList.getColor("power");
        position = propList.getPoint3("position");
    }

    virtual Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample) const override {
        lRec.p = position;
        lRec.wi = (lRec.p - lRec.ref).normalized();
        lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, (lRec.p - lRec.ref).norm() + Epsilon);
        lRec.pdf = 1.f;

        return eval(lRec) / pdf(lRec);
    }

    virtual Color3f eval(const EmitterQueryRecord &lRec) const override {
        return power * INV_FOURPI / (lRec.p - lRec.ref).squaredNorm();
    }

    virtual float pdf(const EmitterQueryRecord &lRec) const override {
        return 1.f;
    }

    virtual std::string toString() const override {
        return tfm::format(
                "PointLight[\n"
                "  power = %s,\n"
                "  position = %s,\n"
                "]",
                power.toString(),
                position.toString());
    }
};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END