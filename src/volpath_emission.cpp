#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/phasefunction.h>
#include <stack>
NORI_NAMESPACE_BEGIN
using namespace std;

class VolPathEmissionIntegrator : public Integrator {
private:
    float ray_length;
public:
    VolPathEmissionIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f color = 0;
        Color3f t = 1.f;
        Ray3f currRay = ray;
        Ray3f trackRay;
        int bounce_cnt = 0;
        const Medium* current_medium = nullptr;
        Intersection its;
        Intersection trackIts;
        stack<Medium*> media;
        float w_mat = 1.f;
        bool has_intersection = scene->rayIntersect(currRay, its);
        while (true) {
            MediumQueryRecord mRec(currRay.o, -currRay.d, its.t);
            if (current_medium && current_medium->sample_freepath(mRec, sampler)) {
                // scattered inside the medium
                PhaseFunctionQueryRecord pRec(mRec.wi);
                current_medium->getPhaseFunction()->sample(pRec, sampler->next2D());   //sample direction to next interaction
                //sample direction to next sampling
                currRay = Ray3f(mRec.p, pRec.wo);
                has_intersection = scene->rayIntersect(currRay, its);
                t *= mRec.ret;
                if (current_medium->isEmissive()) {
                    color += t * mRec.radiance;           
                }
            }
            else {
                if (current_medium) {
                    t *= mRec.ret;
                }
                if (!has_intersection) {
                    break;
                }

                /* compute Le when intersection is an emitter */
                if (its.mesh->isEmitter()) {
                    EmitterQueryRecord lRec(currRay.o, its.p, its.shFrame.n);
                    Color3f Li = its.mesh->getEmitter()->eval(lRec);
                    color += w_mat * t * its.mesh->getEmitter()->eval(lRec);
                }

                //determine whether to sample it or not
                if (scene->sampleEmitter(sampler->next1D())) {
                    /* sample emitter */
                    auto light = scene->getRandomEmitter(sampler->next1D());
                    int n_lights = scene->getLights().size();
                    EmitterQueryRecord lRec(its.p);
                    Color3f Li = light->sample(lRec, sampler->next2D());

                    /* compute w_em */
                    if (!scene->rayIntersect(lRec.shadowRay)) {  // if shadow ray is NOT occluded
                        float cosTheta = its.shFrame.n.dot(lRec.wi);
                        BSDFQueryRecord bRec(its.toLocal(-currRay.d), its.toLocal(lRec.wi), ESolidAngle);
                        bRec.uv = its.uv;
                        float pdf_em = light->pdf(lRec);
                        float pdf_mat = its.mesh->getBSDF()->pdf(bRec);
                        if (pdf_em + pdf_mat > 1e-8) {
                            float w_em = pdf_em / (pdf_mat + pdf_em);
                            color += w_em * t * its.mesh->getBSDF()->eval(bRec) * Li * cosTheta * n_lights;
                        }
                    }
                }
                else {
                    /* sample emissive media*/
                    auto emissive = scene->getRandomEmissiveMedia(sampler->next1D());
                    int n_emissive = scene->getEmissiveMedia().size();
                    MediumQueryRecord ems_mRec(its.p);
                    Color3f Li = emissive->sample_radiance(ems_mRec, sampler);
                    Color3f tr = emissive->Tr(ems_mRec, sampler);
                    Li *= tr;
                    if (!scene->rayIntersect(ems_mRec.shadowRay)) {  // if shadow ray is NOT occluded
                        float cosTheta = its.shFrame.n.dot(ems_mRec.wi);
                        BSDFQueryRecord bRec(its.toLocal(-currRay.d), its.toLocal(ems_mRec.wi), ESolidAngle);
                        bRec.uv = its.uv;
                        float pdf_em = ems_mRec.radiance_pdf;
                        float pdf_mat = its.mesh->getBSDF()->pdf(bRec);
                        if (pdf_em + pdf_mat > 1e-8) {
                            float w_em = pdf_em / (pdf_mat + pdf_em);
                            color += w_em * t * its.mesh->getBSDF()->eval(bRec) * Li * cosTheta * n_emissive;
                        }
                    }
                }

                /* Russian roulette */
                float p = std::min(t.maxCoeff(), 0.99f);
                if (sampler->next1D() > p || p <= 0.f) {
                    break;
                }
                t /= p;

                /* sample brdf */
                BSDFQueryRecord bRec(its.toLocal((-currRay.d).normalized()));
                bRec.p = its.p;
                bRec.uv = its.uv;
                Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

                /* recursion */
                t *= brdf;
                currRay = Ray3f(bRec.p, its.toWorld(bRec.wo));

                /* current pdf mat */
                float pdf_mat = its.mesh->getBSDF()->pdf(bRec);
                if (Frame::cosTheta(its.shFrame.toLocal(currRay.d)) >= 0.0f && its.mesh->getMedium() == current_medium) {
                    current_medium = nullptr;
                }
                else if (Frame::cosTheta(its.shFrame.toLocal(currRay.d)) < 0.0f && its.mesh->getMedium() != nullptr) {
                    current_medium = its.mesh->getMedium();
                }

                has_intersection = scene->rayIntersect(currRay, its);

                /* compute w_mat */
                if (its.mesh->isEmitter()) {
                    EmitterQueryRecord lRec(currRay.o, its.p, its.shFrame.n);
                    float pdf_em = its.mesh->getEmitter()->pdf(lRec);
                    if (pdf_em + pdf_mat > 1e-8) {
                        w_mat = pdf_mat / (pdf_mat + pdf_em);
                    }
                    else {
                        w_mat = 0.f;
                    }
                }

                if (bRec.measure == EDiscrete)
                    w_mat = 1.f;
            }
        }

        return color;
    }

    std::string toString() const {
        return "VolPathEmissiveIntegrator[]";
    }
};

NORI_REGISTER_CLASS(VolPathEmissionIntegrator, "vol_path_ems");
NORI_NAMESPACE_END