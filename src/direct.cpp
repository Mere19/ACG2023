#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class Direct : public Integrator {
public:
    Direct(const PropertyList &props) {

    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* return black when no intersection */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        /* compute total light intensity and intersection point */
        Color3f light_total(0.f);
        auto bsdf = its.mesh->getBSDF();
        auto lights = scene->getLights();
        for (auto const l: lights)
        {
            EmitterQueryRecord lRec(its.p);
            Color3f color = l->sample(lRec, Point2f(0.f,0.f));
            Intersection shIts; // shadow ray intersection
            if(!scene->rayIntersect(lRec.shadowRay, shIts))
            {
                BSDFQueryRecord bRec(its.toLocal(lRec.wi), its.toLocal(-ray.d), ESolidAngle);
                bRec.uv = its.uv;
                light_total += bsdf->eval(bRec) * color * lRec.wi.dot(its.shFrame.n);
            }
        }
        
        return light_total;
    }

    std::string toString() const {
        return "Direct[]";
    }
};

NORI_REGISTER_CLASS(Direct, "direct");
NORI_NAMESPACE_END