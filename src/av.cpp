#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AverageVisibility : public Integrator {
public:
    float length;
    AverageVisibility (const PropertyList &props) {
        this->length = props.getFloat("length");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* find the surface intersected by the camera ray */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(1.0f);

        /* generate a point on the hemisphere and trace a ray into this direction */
        Vector3f n =  Warp::sampleUniformHemisphere(sampler, its.shFrame.n);
        Ray3f r(its.p, n, Epsilon, this->length + Epsilon);
        if (!scene->rayIntersect(r, its)) {
            return Color3f(1.0f);
        } else {
            return Color3f(0.0f);
        }
    }

    std::string toString() const {
        return "AverageVisibility[]";
    }
};

NORI_REGISTER_CLASS(AverageVisibility, "av");
NORI_NAMESPACE_END