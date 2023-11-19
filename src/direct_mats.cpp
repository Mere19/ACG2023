#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class DirectMats : public Integrator {
public:
    DirectMats(const PropertyList &props) {

    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;
        Color3f color = 0;

        /* return black when no intersection */
        if (!scene->rayIntersect(ray, its))
            return 0;

        /* compute Le when intersection is an emitter */
        if (its.mesh->isEmitter()) {
            EmitterQueryRecord lRec(ray.o, its.p, its.shFrame.n);
            color += its.mesh->getEmitter()->eval(lRec);
        }

        /* sample brdf */
        BSDFQueryRecord bRec(its.toLocal((-ray.d).normalized()));
        bRec.p = its.p;
        bRec.uv = its.uv;
        Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

        /* shadow ray pointing to the emitter */
        Ray3f ray_wo(its.p, its.toWorld(bRec.wo), Epsilon, INFINITY);
        Intersection its_wo;
        if (scene->rayIntersect(ray_wo, its_wo)) {
            if (its_wo.mesh->isEmitter()) {
                EmitterQueryRecord lRec(its.p, its_wo.p, its_wo.shFrame.n);
                color += brdf * its_wo.mesh->getEmitter()->eval(lRec);
            }
        }

        return color;
    }

    std::string toString() const {
        return "DirectMats[]";
    }
};

NORI_REGISTER_CLASS(DirectMats, "direct_mats")
NORI_NAMESPACE_END