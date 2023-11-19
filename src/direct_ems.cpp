#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class DirectEms : public Integrator {
public:
    DirectEms(const PropertyList &props) {

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

        /* compute Lr */
        int n_lights = scene->getLights().size();
        auto l = scene->getRandomEmitter(sampler->next1D());
        EmitterQueryRecord lRec(its.p);
        Color3f Li = l->sample(lRec, sampler->next2D());
        if (!scene->rayIntersect(lRec.shadowRay)) {  // if shadow ray is NOT occluded
            float cosTheta = its.shFrame.n.dot(lRec.wi);
            BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(lRec.wi), ESolidAngle);
            color += its.mesh->getBSDF()->eval(bRec) * Li * cosTheta * n_lights;
        }

        return color;
    }

    std::string toString() const {
        return "DirectEms[]";
    }
};

NORI_REGISTER_CLASS(DirectEms, "direct_ems")
NORI_NAMESPACE_END