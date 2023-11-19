#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class DirectMis : public Integrator {
public:
    DirectMis(const PropertyList &props) {

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

        /* initialization */
        auto bsdf = its.mesh->getBSDF();
        auto light = scene->getRandomEmitter(sampler->next1D());
        int n_lights = scene->getLights().size();

        /* sample brdf */
        BSDFQueryRecord bRec(its.toLocal((-ray.d).normalized()));
        bRec.p = its.p;
        bRec.uv = its.uv;
        bRec.measure = ESolidAngle;
        Color3f brdf = bsdf->sample(bRec, sampler->next2D());

        /* compute brdf component */
        Ray3f ray_wo(its.p, its.toWorld(bRec.wo), Epsilon, INFINITY);
        Intersection its_wo;
        if (scene->rayIntersect(ray_wo, its_wo)) {
            if (its_wo.mesh->isEmitter()) {
                EmitterQueryRecord lRec(its.p, its_wo.p, its_wo.shFrame.n);
                Color3f Li = its_wo.mesh->getEmitter()->eval(lRec);
                float pdf_mat = bsdf->pdf(bRec);
                float pdf_em = its_wo.mesh->getEmitter()->pdf(lRec);
                if (pdf_em + pdf_mat > 1e-8) {
                    float w_mat = pdf_mat / (pdf_mat + pdf_em);
                    color += w_mat * brdf * Li;
                }
            }
        }

        /* sample emitter */
        EmitterQueryRecord lRec(its.p);
        Color3f Li = light->sample(lRec, sampler->next2D());

        /* compute emitter component */
        if (!scene->rayIntersect(lRec.shadowRay)) {  // if shadow ray is NOT occluded
            float cosTheta = its.shFrame.n.dot(lRec.wi);
            BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(lRec.wi), ESolidAngle);
            bRec.uv = its.uv;
            float pdf_em = light->pdf(lRec);
            float pdf_mat = bsdf->pdf(bRec);
            if (pdf_em + pdf_mat > 1e-8) {
                float w_em = pdf_em / (pdf_mat + pdf_em);
                color += w_em * bsdf->eval(bRec) * Li * cosTheta * n_lights;
            }
        }

        return color;
    }

    std::string toString() const {
        return "DirectMis[]";
    }
};

NORI_REGISTER_CLASS(DirectMis, "direct_mis")
NORI_NAMESPACE_END