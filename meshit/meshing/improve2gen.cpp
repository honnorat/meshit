#include "improve2.hpp"

#include "global.hpp"

namespace meshit {

    class ImprovementRule
    {
     public:
        Array<Element2d> oldels;
        Array<Element2d> newels;
        Array<INDEX_2> deledges;
        Array<int> incelsonnode;
        Array<int> reused;
        int bonus;
        int onp;
    };

    void MeshOptimize2d::GenericImprove(Mesh& mesh)
    {
        if (!faceindex) {
            if (writestatus)
                MESHIT_LOG_DEBUG("Generic Improve");

            for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
                GenericImprove(mesh);

            faceindex = 0;
        }

        int np = mesh.GetNP();
        int ne = mesh.GetNSE();

        bool ok;
        int olddef, newdef;

        Array<ImprovementRule*> rules;
        Array<size_t> elmap;
        Array<size_t> elrot;
        Array<PointIndex> pmap;
        Array<PointGeomInfo> pgi;

        int surfnr = mesh.GetFaceDescriptor(faceindex).SurfNr();

        ImprovementRule* r1;

        // 2 triangles to quad
        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3));
        r1->oldels.push_back(Element2d(3, 2, 4));
        r1->newels.push_back(Element2d(1, 2, 4, 3));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->onp = 4;
        r1->bonus = 2;
        rules.push_back(r1);

        // 2 quad to 1 quad
        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(4, 3, 2, 5));
        r1->newels.push_back(Element2d(1, 2, 5, 4));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->deledges.push_back(INDEX_2(3, 4));
        r1->onp = 5;
        r1->bonus = 0;
        rules.push_back(r1);

        // std::swap quads
        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(3, 2, 5, 6));
        r1->newels.push_back(Element2d(1, 6, 3, 4));
        r1->newels.push_back(Element2d(1, 2, 5, 6));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->onp = 6;
        r1->bonus = 0;
        rules.push_back(r1);

        // three quads to 2
        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(2, 5, 6, 3));
        r1->oldels.push_back(Element2d(3, 6, 7, 4));
        r1->newels.push_back(Element2d(1, 2, 5, 4));
        r1->newels.push_back(Element2d(4, 5, 6, 7));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->deledges.push_back(INDEX_2(3, 4));
        r1->deledges.push_back(INDEX_2(3, 6));
        r1->onp = 7;
        r1->bonus = -1;
        rules.push_back(r1);

        // quad + 2 connected trigs to quad
        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(2, 5, 3));
        r1->oldels.push_back(Element2d(3, 5, 4));
        r1->newels.push_back(Element2d(1, 2, 5, 4));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->deledges.push_back(INDEX_2(3, 4));
        r1->deledges.push_back(INDEX_2(3, 5));
        r1->onp = 5;
        r1->bonus = 0;
        rules.push_back(r1);

        // quad + 2 non-connected trigs to quad (a and b)
        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(2, 6, 3));
        r1->oldels.push_back(Element2d(1, 4, 5));
        r1->newels.push_back(Element2d(1, 3, 4, 5));
        r1->newels.push_back(Element2d(1, 2, 6, 3));
        r1->deledges.push_back(INDEX_2(1, 4));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->onp = 6;
        r1->bonus = 0;
        rules.push_back(r1);

        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(2, 6, 3));
        r1->oldels.push_back(Element2d(1, 4, 5));
        r1->newels.push_back(Element2d(1, 2, 4, 5));
        r1->newels.push_back(Element2d(4, 2, 6, 3));
        r1->deledges.push_back(INDEX_2(1, 4));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->onp = 6;
        r1->bonus = 0;
        rules.push_back(r1);

        // two quad + trig -> one quad + trig
        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(2, 5, 6, 3));
        r1->oldels.push_back(Element2d(4, 3, 6));
        r1->newels.push_back(Element2d(1, 2, 6, 4));
        r1->newels.push_back(Element2d(2, 5, 6));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->deledges.push_back(INDEX_2(3, 4));
        r1->deledges.push_back(INDEX_2(3, 6));
        r1->onp = 6;
        r1->bonus = -1;
        rules.push_back(r1);

        // std::swap quad + trig (a and b)
        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(2, 5, 3));
        r1->newels.push_back(Element2d(2, 5, 3, 4));
        r1->newels.push_back(Element2d(1, 2, 4));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->onp = 5;
        r1->bonus = 0;
        rules.push_back(r1);

        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(2, 5, 3));
        r1->newels.push_back(Element2d(1, 2, 5, 3));
        r1->newels.push_back(Element2d(1, 3, 4));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->onp = 5;
        r1->bonus = 0;
        rules.push_back(r1);


        // 2 quads to quad + 2 trigs
        r1 = new ImprovementRule;
        r1->oldels.push_back(Element2d(1, 2, 3, 4));
        r1->oldels.push_back(Element2d(3, 2, 5, 6));
        r1->newels.push_back(Element2d(1, 5, 6, 4));
        r1->newels.push_back(Element2d(1, 2, 5));
        r1->newels.push_back(Element2d(4, 6, 3));
        r1->deledges.push_back(INDEX_2(2, 3));
        r1->onp = 6;
        r1->bonus = 0;
        //    rules.Append (r1);

        Array<int> mapped(rules.size());
        Array<int> used(rules.size());
        used = 0;
        mapped = 0;

        for (size_t ri = 0; ri < rules.size(); ri++) {
            ImprovementRule& rule = *rules[ri];
            rule.incelsonnode.resize(rule.onp);
            rule.reused.resize(rule.onp);

            for (int j = 0; j < rule.onp; j++) {
                rule.incelsonnode[j] = 0;
                rule.reused[j] = 0;
            }

            for (size_t j = 0; j < rule.oldels.size(); j++) {
                const Element2d& el = rule.oldels[j];
                for (size_t k = 1; k <= el.GetNP(); k++) {
                    rule.incelsonnode[el.PNum(k) - 1]--;
                }
            }

            for (size_t j = 0; j < rule.newels.size(); j++) {
                const Element2d& el = rule.newels[j];
                for (size_t k = 1; k <= el.GetNP(); k++) {
                    rule.incelsonnode[el.PNum(k) - 1]++;
                    rule.reused[el.PNum(k) - 1] = 1;
                }
            }
        }

        TABLE<int> elonnode(np);
        Array<int> nelonnode(np);
        TABLE<SurfaceElementIndex> nbels(ne);

        nelonnode = -4;

        for (SurfaceElementIndex sei = 0; sei < ne; sei++) {
            const Element2d& el = mesh.SurfaceElement(sei);

            if (el.GetIndex() == faceindex && !el.IsDeleted()) {
                for (size_t j = 0; j < el.GetNP(); j++)
                    elonnode.Add(el[j], sei);
            }
            if (!el.IsDeleted()) {
                for (size_t j = 0; j < el.GetNP(); j++)
                    nelonnode[el[j]]++;
            }
        }

        for (SurfaceElementIndex sei = 0; sei < ne; sei++) {
            const Element2d& el = mesh.SurfaceElement(sei);
            if (el.GetIndex() == faceindex && !el.IsDeleted()) {
                for (size_t j = 0; j < el.GetNP(); j++) {
                    for (size_t k = 0; k < elonnode[el[j]].size(); k++) {
                        int nbel = elonnode[el[j]][k];
                        bool inuse = false;
                        for (size_t l = 0; l < nbels[sei].size(); l++) {
                            if (nbels[sei][l] == nbel) {
                                inuse = true;
                            }
                        }
                        if (!inuse)
                            nbels.Add(sei, nbel);
                    }
                }
            }
        }
        for (size_t ri = 0; ri < rules.size(); ri++) {
            const ImprovementRule& rule = *rules[ri];

            elmap.resize(rule.oldels.size());
            elrot.resize(rule.oldels.size());
            pmap.resize(rule.onp);
            pgi.resize(rule.onp);

            for (size_t sei = 0; sei < ne; sei++) {

                if (mesh.SurfaceElement(sei+1).IsDeleted()) continue;

                elmap[0] = sei;
                FlatArray<SurfaceElementIndex> neighbours = nbels[sei];

                for (elrot[0] = 0; elrot[0] < mesh.SurfaceElement(sei+1).GetNP(); elrot[0]++) {
                    const Element2d& el0 = mesh.SurfaceElement(sei+1);
                    const Element2d& rel0 = rule.oldels[0];

                    if (el0.GetIndex() != faceindex) continue;
                    if (el0.IsDeleted()) continue;
                    if (el0.GetNP() != rel0.GetNP()) continue;


                    pmap = PointIndex(-1);

                    for (size_t k = 0; k < el0.GetNP(); k++) {
                        pmap[rel0[k] - 1] = el0.PNumMod(k + elrot[0] + 1);
                        pgi[rel0[k] - 1] = el0.GeomInfoPiMod(k + elrot[0] + 1);
                    }

                    ok = 1;
                    for (size_t i = 1; i < elmap.size(); i++) {
                        // try to find a mapping for reference-element i
                        const Element2d& rel = rule.oldels[i];
                        bool possible = 0;

                        for (elmap[i] = 0; elmap[i] < neighbours.size(); elmap[i]++) {
                            const Element2d& el = mesh.SurfaceElement(neighbours[elmap[i]]);
                            if (el.IsDeleted()) continue;
                            if (el.GetNP() != rel.GetNP()) continue;

                            for (elrot[i] = 0; elrot[i] < rel.GetNP(); elrot[i]++) {
                                possible = 1;
                                for (size_t k = 0; k < rel.GetNP(); k++) {
                                    if (pmap[rel[k] - 1] != -1 && pmap[rel[k] - 1] != el.PNumMod(k + elrot[i] + 1)) {
                                        possible = 0;
                                    }
                                }
                                if (possible) {
                                    for (size_t k = 0; k < el.GetNP(); k++) {
                                        pmap[rel[k] - 1] = el.PNumMod(k + elrot[i] + 1);
                                        pgi[rel[k] - 1] = el.GeomInfoPiMod(k + elrot[i] + 1);
                                    }
                                    break;
                                }
                            }
                            if (possible) break;
                        }

                        if (!possible) {
                            ok = 0;
                            break;
                        }

                        elmap[i] = neighbours[elmap[i]];
                    }

                    for (size_t i = 0; ok && i < rule.deledges.size(); i++) {
                        ok = !mesh.IsSegment(pmap[rule.deledges[i].I1() - 1],
                                             pmap[rule.deledges[i].I2() - 1]);
                    }

                    if (!ok) continue;

                    mapped[ri]++;

                    olddef = 0;
                    for (size_t j = 0; j < pmap.size(); j++) {
                        int ii = nelonnode[pmap[j]];
                        olddef += ii * ii;
                    }
                    olddef += rule.bonus;

                    newdef = 0;
                    for (size_t j = 0; j < pmap.size(); j++) {
                        if (rule.reused[j]) {
                            int ii = nelonnode[pmap[j]] + rule.incelsonnode[j];
                            newdef += ii * ii;
                        }
                    }

                    if (newdef > olddef)
                        continue;

                    // calc metric badness
                    double bad1 = 0, bad2 = 0;
                    Vec3d n;

                    GetNormalVector(surfnr, mesh.Point(pmap[0]), pgi[0], n);

                    for (size_t j = 0; j < rule.oldels.size(); j++) {
                        bad1 += mesh.SurfaceElement(elmap[j]+1).CalcJacobianBadness(mesh.Points(), n);
                    }
                    // check new element:
                    for (size_t j = 0; j < rule.newels.size(); j++) {
                        const Element2d& rnel = rule.newels[j];
                        Element2d nel(rnel.GetNP());
                        for (size_t k = 1; k <= rnel.GetNP(); k++) {
                            nel.PNum(k) = pmap[rnel.PNum(k) - 1];
                        }
                        bad2 += nel.CalcJacobianBadness(mesh.Points(), n);
                    }

                    if (bad2 > 1e3) continue;

                    if (newdef == olddef && bad2 > bad1) continue;

                    // generate new element:
                    for (size_t j = 0; j < rule.newels.size(); j++) {
                        const Element2d& rnel = rule.newels[j];
                        Element2d nel(rnel.GetNP());
                        nel.SetIndex(faceindex);
                        for (size_t k = 1; k <= rnel.GetNP(); k++) {
                            nel.PNum(k) = pmap[rnel.PNum(k) - 1];
                            nel.GeomInfoPi(k) = pgi[rnel.PNum(k) - 1];
                        }
                        mesh.AddSurfaceElement(nel);
                    }

                    for (size_t j = 0; j < rule.oldels.size(); j++)
                        mesh.DeleteSurfaceElement(elmap[j]);

                    for (size_t j = 0; j < pmap.size(); j++)
                        nelonnode[pmap[j]] += rule.incelsonnode[j];

                    used[ri]++;
                }
            }
        }

        mesh.Compress();

        for (size_t ri = 0; ri < rules.size(); ri++) {
            MESHIT_LOG_DEBUG("rule " << ri + 1 << " " << mapped[ri] << "/" << used[ri] << " mapped/used");
        }
    }
}  // namespace meshit
