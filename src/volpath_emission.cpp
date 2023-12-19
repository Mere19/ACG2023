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

                //emission sampling
                auto emissive = scene->getRandomEmissiveMedia(sampler->next1D());
                int n_emissive = scene->getEmissiveMedia().size();
                MediumQueryRecord mRec_ems(its.p);
                Color3f Li = emissive->sample_radiance(mRec_ems, sampler);
                float pdf_em = Epsilon;
                Ray3f testRay = Ray3f(mRec_ems.p, (emissive->getShape()->getBoundingBox().getCenter() - mRec_ems.p).normalized(), Epsilon, (emissive->getShape()->getBoundingBox().getCenter() - mRec_ems.p).norm());
                Ray3f sampleRay = Ray3f(its.p, (mRec_ems.p - its.p).normalized(), Epsilon, (mRec_ems.p - its.p).norm());

                if (!scene->rayCurrIntersect(testRay, emissive->getShape())) {
                    //if not rejected (inside bbox, outside the medium), accept the distribution
                    if (!scene->rayIntersect(sampleRay)) {
                        //if not occuluded set the pdf, otherwise 0 contribution.
                        pdf_em = mRec_ems.radiance_pdf;
                    }
                }

                Color3f Transmittance = transmittance(scene, sampler, sampleRay, current_medium);
                PhaseFunctionQueryRecord pRec_ems(mRec.wi, sampleRay.d, ESolidAngle);
                float pdf_mat = current_medium->getPhaseFunction()->pdf(pRec_ems);
                if (pdf_em >= 2 * Epsilon) {
                    //ignore invalid contributions
                    color += t * Li * n_emissive / pdf_em;
                }


                // scattered inside the medium
                PhaseFunctionQueryRecord pRec(mRec.wi);
                current_medium->getPhaseFunction()->sample(pRec, sampler->next2D());   //sample direction to next interaction
                //sample direction to next sampling
                currRay = Ray3f(mRec.p, pRec.wo);
                has_intersection = scene->rayIntersect(currRay, its);
                t *= mRec.ret;
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

    /* assumption: trRay.o, trRay.d, trRay.maxt*/
    Color3f transmittance(const Scene* scene, Sampler* sampler, const Ray3f & trRay, const Medium* medium) const {
        float start = 0.f;
        Color3f Tr = Color3f(1);
        const Medium* current_medium = medium;
        Ray3f currRay = trRay;
        Intersection its;
        int count = 0;
        int maxDepth = 30;
        //based on the assumption that no overlapping medium
        while (true) {
            count++;
            if (count > maxDepth) {
                break;
            }
            if (start >= trRay.maxt) {
                break;
            }
            bool intersect = scene->rayIntersect(currRay, its);
            if (!intersect|| its.t >= currRay.maxt)  {
                if (current_medium) {
                    MediumQueryRecord segMRec(currRay.o, currRay(currRay.maxt));
                    return Tr * current_medium->Tr(segMRec, sampler);
                }
                else {
                    return Tr;
                }
            }
            else {
                if (!its.mesh->isMedium()) {
                    return Color3f(0.0f);
                }
                else {
                    if (Frame::cosTheta(its.shFrame.toLocal(currRay.d)) >= 0.0f && current_medium) {
                        MediumQueryRecord segMRec(currRay.o, currRay(currRay.maxt));
                        Tr *= current_medium->Tr(segMRec, sampler);
                        current_medium = nullptr;
                        start += its.t;
                        currRay.maxt = trRay.maxt - start;
                        currRay.o = trRay(start);
                        currRay.update();
                    }
                    if (Frame::cosTheta(its.shFrame.toLocal(currRay.d)) < 0.0f && !current_medium) {
                        current_medium = its.mesh->getMedium();
                        start += its.t;
                        currRay.maxt = trRay.maxt - start;
                        currRay.o = trRay(start);
                        currRay.update();
                    }
                }
            }

        }
    }

    std::string toString() const {
        return "VolPathEmissiveIntegrator[]";
    }
};

NORI_REGISTER_CLASS(VolPathEmissionIntegrator, "vol_path_ems");
NORI_NAMESPACE_END