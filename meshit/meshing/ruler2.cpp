#include "meshing2.hpp"

namespace meshit
{
    static double CalcElementBadness(const std::vector<Point2d>& points, const Element2d& elem)
    {
        // badness = sqrt(3) /36 * circumference^2 / area - 1 +
        //           h / li + li / h - 2

        Vec2d v12, v13, v23;
        double l12, l13, l23, cir, area;
        static const double c = sqrt(3.0) / 36;

        v12 = points[elem.PointID(1) - 1] - points[elem.PointID(0) - 1];
        v13 = points[elem.PointID(2) - 1] - points[elem.PointID(0) - 1];
        v23 = points[elem.PointID(2) - 1] - points[elem.PointID(1) - 1];

        l12 = v12.Length();
        l13 = v13.Length();
        l23 = v23.Length();

        cir = l12 + l13 + l23;
        area = 0.5 * (v12.X() * v13.Y() - v12.Y() * v13.X());
        if (area < 1e-6) {
            return 1e8;
        }

        return 10 * (c * cir * cir / area - 1)
               + 1 / l12 + l12 + 1 / l13 + l13 + 1 / l23 + l23 - 6;
    }

    int Meshing2::ApplyRules(std::vector<Point2d>& lpoints,
                             std::vector<int>& legalpoints,
                             size_t maxlegalpoint,
                             std::vector<INDEX_2>& llines1,
                             size_t maxlegalline,
                             std::vector<Element2d>& elements,
                             std::vector<uint32_t>& dellines,
                             int tolerance,
                             const MeshingParameters& mp)
    {
        double maxerr = 0.5 + 0.3 * tolerance;
        double minelerr = 2 + 0.5 * tolerance * tolerance;

        size_t noldlp = lpoints.size();
        size_t noldll = llines1.size();

        std::vector<Point2d> tempnewpoints;
        std::vector<INDEX_2> tempnewlines;
        std::vector<uint32_t> tempdellines;
        std::vector<Element2d> tempelements;

        elements.resize(0);
        dellines.resize(0);

        // check every rule

        size_t found = 0;  // rule number

        std::vector<uint32_t> pnearness(noldlp, 1000);
        std::vector<uint32_t> lnearness(llines1.size());

        pnearness[llines1[0][0] - 1] = 0;
        pnearness[llines1[0][1] - 1] = 0;

        const size_t MAX_NEARNESS = 3;

        for (size_t cnt = 0; cnt < MAX_NEARNESS; cnt++) {
            bool ok = true;
            for (size_t i = 0; i < maxlegalline; i++) {
                const INDEX_2& hline = llines1[i];

                size_t minn = std::min(pnearness[hline[0] - 1],
                                       pnearness[hline[1] - 1]);

                for (size_t j = 0; j < 2; j++) {
                    if (pnearness[hline[j] - 1] > minn + 1) {
                        ok = false;
                        pnearness[hline[j] - 1] = minn + 1;
                    }
                }
            }
            if (!ok) break;
        }

        for (size_t i = 0; i < maxlegalline; i++) {
            lnearness[i] = pnearness[llines1[i][0] - 1] + pnearness[llines1[i][1] - 1];
        }

        // resort lines after lnearness
        std::vector<INDEX_2> llines(llines1.size());
        std::vector<size_t> sortlines(llines1.size());
        uint32_t lnearness_class[MAX_NEARNESS];

        for (size_t j = 0; j < MAX_NEARNESS; j++) {
            lnearness_class[j] = 0;
        }
        for (size_t i = 0; i < maxlegalline; i++) {
            if (lnearness[i] < MAX_NEARNESS)
                lnearness_class[lnearness[i]]++;
        }
        int cumm = 0;
        for (size_t j = 0; j < MAX_NEARNESS; j++) {
            int hcnt = lnearness_class[j];
            lnearness_class[j] = cumm;
            cumm += hcnt;
        }

        for (size_t i = 0; i < maxlegalline; i++)
            if (lnearness[i] < MAX_NEARNESS) {
                llines[lnearness_class[lnearness[i]]] = llines1[i];
                sortlines[lnearness_class[lnearness[i]]] = i + 1;
                lnearness_class[lnearness[i]]++;
            } else {
                llines[cumm] = llines1[i];
                sortlines[cumm] = i + 1;
                cumm++;
            }

        for (size_t i = maxlegalline; i < llines1.size(); i++) {
            llines[cumm] = llines1[i];
            sortlines[cumm] = i + 1;
            cumm++;
        }

        for (size_t i = 0; i < maxlegalline; i++) {
            lnearness[i] = pnearness[llines[i][0] - 1] + pnearness[llines[i][1] - 1];
        }

        std::vector<int> pused(maxlegalpoint, 0);
        std::vector<bool> lused(maxlegalline, false);
        std::vector<int> pmap;
        std::vector<bool> pfixed;
        std::vector<size_t> lmap;

        for (size_t ri = 0; ri < rules.size(); ri++) {
            netrule* rule = rules[ri];

            if (rule->GetQuality() > tolerance) continue;

            pmap.assign(rule->GetNP(), 0);
            lmap.assign(rule->GetNL(), 0);

            lused[0] = true;
            lmap[0] = 1;

            for (int j = 0; j < 2; j++) {
                pmap[rule->GetLine(0)[j] - 1] = llines[0][j];
                pused[llines[0][j] - 1]++;
            }

            size_t nlok = 1;
            bool ok = false;

            while (nlok > 0) {
                if (nlok < rule->GetNOldL()) {
                    ok = 0;

                    size_t maxline = (rule->GetLNearness(nlok) < MAX_NEARNESS)
                                     ? lnearness_class[rule->GetLNearness(nlok)] : maxlegalline;

                    while (!ok && lmap[nlok] < maxline) {
                        lmap[nlok]++;
                        int locli = lmap[nlok];

                        if (lnearness[locli - 1] > rule->GetLNearness(nlok)) continue;
                        if (lused[locli - 1]) continue;

                        ok = 1;

                        INDEX_2 loclin = llines[locli - 1];
                        Vec2d linevec = lpoints[loclin.I2() - 1] - lpoints[loclin.I1() - 1];

                        if (rule->CalcLineError(nlok, linevec) > maxerr) {
                            ok = 0;
                            continue;
                        }

                        for (int j = 0; j < 2; j++) {
                            int refpi = rule->GetLine(nlok)[j] - 1;

                            if (pmap[refpi] != 0) {
                                if (pmap[refpi] != loclin[j]) {
                                    ok = 0;
                                    break;
                                }
                            } else {
                                if (rule->CalcPointDist(refpi, lpoints[loclin[j] - 1]) > maxerr
                                    || !legalpoints[loclin[j] - 1]
                                    || pused[loclin[j] - 1]) {
                                    ok = 0;
                                    break;
                                }
                            }
                        }
                    }

                    if (ok) {
                        int locli = lmap[nlok];
                        INDEX_2 loclin = llines[locli - 1];

                        lused[locli - 1] = true;
                        for (int j = 0; j < 2; j++) {
                            pmap[rule->GetLine(nlok)[j] - 1] = loclin[j];
                            pused[loclin[j] - 1]++;
                        }
                        nlok++;
                    } else {
                        lmap[nlok] = 0;
                        nlok--;

                        lused[lmap[nlok] - 1] = false;
                        for (int j = 0; j < 2; j++) {
                            pused[llines[lmap[nlok] - 1][j] - 1]--;
                            if (!pused[llines[lmap[nlok] - 1][j] - 1]) {
                                pmap[rule->GetLine(nlok)[j] - 1] = 0;
                            }
                        }
                    }
                } else {
                    // all lines are mapped !!
                    // map also all points:

                    int npok = 0;
                    int incnpok = 1;

                    pfixed.resize(pmap.size());
                    for (size_t i = 0; i < pmap.size(); i++) {
                        pfixed[i] = (pmap[i] >= 1);
                    }

                    while (npok >= 0) {
                        if (npok < static_cast<int>(rule->GetNOldP())) {
                            if (pfixed[npok]) {
                                if (incnpok)
                                    npok++;
                                else
                                    npok--;
                            } else {
                                ok = 0;

                                if (pmap[npok]) {
                                    pused[pmap[npok] - 1]--;
                                }
                                while (!ok && pmap[npok] < static_cast<int>(maxlegalpoint)) {
                                    ok = 1;

                                    pmap[npok]++;

                                    if (pused[pmap[npok] - 1]) {
                                        ok = 0;
                                    } else {
                                        if (rule->CalcPointDist(npok, lpoints[pmap[npok] - 1]) > maxerr
                                            || !legalpoints[pmap[npok] - 1]) {
                                            ok = 0;
                                        }
                                    }
                                }

                                if (ok) {
                                    pused[pmap[npok] - 1]++;
                                    npok++;
                                    incnpok = 1;
                                } else {
                                    pmap[npok] = 0;
                                    npok--;
                                    incnpok = 0;
                                }
                            }
                        } else {
                            npok = rule->GetNOldP() - 1;
                            incnpok = 0;

                            if (ok)
                                foundmap[ri]++;

                            ok = 1;

                            // check orientations
                            for (size_t i = 0; i < rule->GetNOrientations(); i++) {
                                if (CW(lpoints[pmap[rule->GetOrientation(i).i1 - 1] - 1],
                                       lpoints[pmap[rule->GetOrientation(i).i2 - 1] - 1],
                                       lpoints[pmap[rule->GetOrientation(i).i3 - 1] - 1])) {
                                    ok = 0;
                                    break;
                                }
                            }

                            if (!ok) continue;

                            Vector oldu(2 * rule->GetNOldP());

                            for (size_t i = 0; i < rule->GetNOldP(); i++) {
                                Vec2d ui(rule->GetPoint(i), lpoints[pmap[i] - 1]);
                                oldu(2 * i + 0) = ui.X();
                                oldu(2 * i + 1) = ui.Y();
                            }

                            rule->SetFreeZoneTransformation(oldu, tolerance);

                            if (!rule->ConvexFreeZone()) {
                                ok = 0;
                            }
                            if (!ok) continue;

                            for (size_t i = 0; i < maxlegalpoint; i++) {
                                if (!pused[i] && rule->IsInFreeZone(lpoints[i])) {
                                    ok = 0;
                                    break;
                                }
                            }
                            if (!ok) continue;

                            for (size_t i = maxlegalpoint; i < lpoints.size(); i++) {
                                if (rule->IsInFreeZone(lpoints[i])) {
                                    ok = 0;
                                    break;
                                }
                            }
                            if (!ok) continue;

                            for (size_t i = 0; i < maxlegalline; i++) {
                                if (!lused[i] && rule->IsLineInFreeZone(lpoints[llines[i].I1() - 1],
                                                                        lpoints[llines[i].I2() - 1])) {
                                    ok = 0;
                                    break;
                                }
                            }
                            if (!ok) continue;

                            for (size_t i = maxlegalline; i < llines.size(); i++) {
                                if (rule->IsLineInFreeZone(lpoints[llines[i].I1() - 1],
                                                           lpoints[llines[i].I2() - 1])) {
                                    ok = 0;
                                    break;
                                }
                            }

                            if (!ok) continue;

                            // Setze neue Punkte:
                            if (rule->GetNOldP() < rule->GetNP()) {
                                Vector newu(rule->GetOldUToNewU().Height());
                                rule->GetOldUToNewU().Mult(oldu, newu);

                                size_t oldnp = rule->GetNOldP();
                                for (size_t i = oldnp; i < rule->GetNP(); i++) {
                                    Point2d np = rule->GetPoint(i);
                                    np.X() += newu(2 * (i - oldnp));
                                    np.Y() += newu(2 * (i - oldnp) + 1);

                                    lpoints.push_back(np);
                                    pmap[i] = lpoints.size();
                                }
                            }

                            // Setze neue Linien:
                            for (size_t i = rule->GetNOldL(); i < rule->GetNL(); i++) {
                                llines.push_back(INDEX_2(pmap[rule->GetLine(i)[0] - 1],
                                                         pmap[rule->GetLine(i)[1] - 1]));
                            }

                            // delete old lines:
                            for (size_t i = 0; i < rule->GetNDelL(); i++) {
                                dellines.push_back(sortlines[lmap[rule->GetDelLine(i) - 1] - 1]);
                            }

                            // insert new elements:
                            for (size_t i = 0; i < rule->GetNE(); i++) {
                                elements.push_back(rule->GetElement(i));
                                for (size_t j = 0; j < 3; j++) {
                                    elements[i].PointID(j) = pmap[elements[i].PointID(j) - 1];
                                }
                            }

                            double elerr = 0;
                            for (size_t i = 0; i < elements.size(); i++) {
                                double hf = CalcElementBadness(lpoints, elements[i]);
                                if (hf > elerr) elerr = hf;
                            }

                            canuse[ri]++;

                            if (elerr < 0.99 * minelerr) {
                                minelerr = elerr;
                                found = ri + 1;

                                tempnewpoints.assign(lpoints.begin() + noldlp, lpoints.end());
                                tempnewlines.assign(llines.begin() + noldll, llines.end());
                                tempdellines = dellines;
                                tempelements = elements;
                            }

                            lpoints.resize(noldlp);
                            llines.resize(noldll);
                            dellines.resize(0);
                            elements.resize(0);
                            ok = 0;
                        }
                    }

                    nlok = rule->GetNOldL() - 1;

                    lused[lmap[nlok] - 1] = false;

                    int refpi1 = rule->GetPointNr1(nlok);
                    int refpi2 = rule->GetPointNr2(nlok);
                    if (--pused[pmap[refpi1] - 1] == 0) pmap[refpi1] = 0;
                    if (--pused[pmap[refpi2] - 1] == 0) pmap[refpi2] = 0;
                }
            }
        }

        if (found) {
            lpoints.insert(lpoints.end(), tempnewpoints.begin(), tempnewpoints.end());
            llines1.insert(llines1.end(), tempnewlines.begin(), tempnewlines.end());
            dellines.insert(dellines.end(), tempdellines.begin(), tempdellines.end());
            elements.insert(elements.end(), tempelements.begin(), tempelements.end());
        }

        return found;
    }
}  // namespace meshit
