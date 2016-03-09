#include "meshing2.hpp"

namespace meshit {

    static double CalcElementBadness(
            const Array<Point2d>& points,
            const Element2d& elem)
    {
        // badness = sqrt(3) /36 * circumference^2 / area - 1 +
        //           h / li + li / h - 2

        Vec2d v12, v13, v23;
        double l12, l13, l23, cir, area;
        static const double c = sqrt(3.0) / 36;

        v12 = points.Get(elem.PNum(2)) - points.Get(elem.PNum(1));
        v13 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(1));
        v23 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(2));

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

    int Meshing2::ApplyRules(
            Array<Point2d>& lpoints,
            Array<int>& legalpoints,
            int maxlegalpoint,
            Array<INDEX_2>& llines1,
            int maxlegalline,
            Array<Element2d>& elements,
            Array<INDEX>& dellines, int tolerance,
            const MeshingParameters& mp)
    {
        double maxerr = 0.5 + 0.3 * tolerance;
        double minelerr = 2 + 0.5 * tolerance * tolerance;

        int noldlp = lpoints.size();
        int noldll = llines1.size();

        ArrayMem<int, 100> pused(maxlegalpoint), lused(maxlegalline);
        ArrayMem<int, 100> pnearness(noldlp), lnearness(llines1.size());
        ArrayMem<int, 20> pmap, pfixed, lmap;

        ArrayMem<Point2d, 100> tempnewpoints;
        ArrayMem<INDEX_2, 100> tempnewlines;
        ArrayMem<int, 100> tempdellines;
        ArrayMem<Element2d, 100> tempelements;

        elements.resize(0);
        dellines.resize(0);

        // check every rule

        int found = 0;  // rule number

        pnearness = 1000;

        for (int j = 0; j < 2; j++) {
            pnearness[llines1[0][j] - 1] = 0;
        }

        const int MAX_NEARNESS = 3;

        for (int cnt = 0; cnt < MAX_NEARNESS; cnt++) {
            bool ok = true;
            for (int i = 0; i < maxlegalline; i++) {
                const INDEX_2& hline = llines1[i];

                int minn = std::min(pnearness.Get(hline[0]), pnearness.Get(hline[1]));

                for (int j = 0; j < 2; j++) {
                    if (pnearness.Get(hline[j]) > minn + 1) {
                        ok = false;
                        pnearness[hline[j]-1] = minn + 1;
                    }
                }
            }
            if (!ok) break;
        }

        for (int i = 0; i < maxlegalline; i++) {
            lnearness[i] = pnearness.Get(llines1[i][0]) + pnearness.Get(llines1[i][1]);
        }

        // resort lines after lnearness
        Array<INDEX_2> llines(llines1.size());
        Array<int> sortlines(llines1.size());
        int lnearness_class[MAX_NEARNESS];

        for (int j = 0; j < MAX_NEARNESS; j++) {
            lnearness_class[j] = 0;
        }
        for (int i = 0; i < maxlegalline; i++) {
            if (lnearness[i] < MAX_NEARNESS)
                lnearness_class[lnearness[i]]++;
        }
        int cumm = 0;
        for (int j = 0; j < MAX_NEARNESS; j++) {
            int hcnt = lnearness_class[j];
            lnearness_class[j] = cumm;
            cumm += hcnt;
        }

        for (int i = 0; i < maxlegalline; i++)
            if (lnearness[i] < MAX_NEARNESS) {
                llines[lnearness_class[lnearness[i]]] = llines1[i];
                sortlines[lnearness_class[lnearness[i]]] = i + 1;
                lnearness_class[lnearness[i]]++;
            } else {
                llines[cumm] = llines1[i];
                sortlines[cumm] = i + 1;
                cumm++;
            }

        for (int i = maxlegalline; i < llines1.size(); i++) {
            llines[cumm] = llines1[i];
            sortlines[cumm] = i + 1;
            cumm++;
        }

        for (int i = 0; i < maxlegalline; i++)
            lnearness[i] = pnearness.Get(llines[i][0]) + pnearness.Get(llines[i][1]);

        lused = 0;
        pused = 0;

        for (int ri = 1; ri <= rules.size(); ri++) {
            netrule* rule = rules.Get(ri);

            if (rule->GetQuality() > tolerance) continue;

            pmap.resize(rule->GetNP());
            lmap.resize(rule->GetNL());

            pmap = 0;
            lmap = 0;

            lused[0] = 1;
            lmap[0] = 1;

            for (int j = 0; j < 2; j++) {
                pmap[rule->GetLine(1)[j] - 1] = llines[0][j];
                pused[llines[0][j] - 1]++;
            }

            int nlok = 2;
            bool ok = false;

            while (nlok >= 2) {
                if (nlok <= static_cast<int>(rule->GetNOldL())) {
                    ok = 0;

                    int maxline = (rule->GetLNearness(nlok) < MAX_NEARNESS)
                                  ? lnearness_class[rule->GetLNearness(nlok)] : maxlegalline;
                    // int maxline = maxlegalline;

                    while (!ok && lmap.Get(nlok) < maxline) {
                        lmap[nlok - 1]++;
                        int locli = lmap.Get(nlok);

                        if (lnearness.Get(locli) > rule->GetLNearness(nlok)) continue;
                        if (lused.Get(locli)) continue;

                        ok = 1;

                        INDEX_2 loclin = llines.Get(locli);
                        Vec2d linevec = lpoints.Get(loclin.I2()) - lpoints.Get(loclin.I1());

                        if (rule->CalcLineError(nlok, linevec) > maxerr) {
                            ok = 0;
                            continue;
                        }

                        for (int j = 0; j < 2; j++) {
                            int refpi = rule->GetLine(nlok)[j];

                            if (pmap.Get(refpi) != 0) {
                                if (pmap.Get(refpi) != loclin[j]) {
                                    ok = 0;
                                    break;
                                }
                            } else {
                                if (rule->CalcPointDist(refpi, lpoints.Get(loclin[j])) > maxerr
                                    || !legalpoints.Get(loclin[j])
                                    || pused.Get(loclin[j])) {
                                    ok = 0;
                                    break;
                                }
                            }
                        }
                    }

                    if (ok) {
                        int locli = lmap.Get(nlok);
                        INDEX_2 loclin = llines.Get(locli);

                        lused[locli - 1] = 1;
                        for (int j = 0; j < 2; j++) {
                            pmap[rule->GetLine(nlok)[j]-1] = loclin[j];
                            pused[loclin[j] - 1]++;
                        }
                        nlok++;
                    } else {
                        lmap[nlok - 1] = 0;
                        nlok--;

                        lused[lmap.Get(nlok) - 1] = 0;
                        for (int j = 0; j < 2; j++) {
                            pused[llines.Get(lmap.Get(nlok))[j] - 1]--;
                            if (!pused.Get(llines.Get(lmap.Get(nlok))[j]))
                                pmap[rule->GetLine(nlok)[j]-1] = 0;
                        }
                    }
                } else {
                    // all lines are mapped !!
                    // map also all points:

                    int npok = 1;
                    int incnpok = 1;

                    pfixed.resize(pmap.size());
                    for (int i = 0; i < pmap.size(); i++) {
                        pfixed[i] = (pmap[i] >= 1);
                    }

                    while (npok >= 1) {
                        if (npok <= static_cast<int>(rule->GetNOldP())) {
                            if (pfixed[npok - 1]) {
                                if (incnpok)
                                    npok++;
                                else
                                    npok--;
                            } else {
                                ok = 0;

                                if (pmap.Get(npok))
                                    pused[pmap.Get(npok) - 1]--;

                                while (!ok && pmap.Get(npok) < maxlegalpoint) {
                                    ok = 1;

                                    pmap[npok - 1]++;

                                    if (pused.Get(pmap.Get(npok))) {
                                        ok = 0;
                                    } else {
                                        if (rule->CalcPointDist(npok, lpoints.Get(pmap.Get(npok))) > maxerr
                                            || !legalpoints.Get(pmap.Get(npok)))

                                            ok = 0;
                                    }
                                }

                                if (ok) {
                                    pused[pmap.Get(npok) - 1]++;
                                    npok++;
                                    incnpok = 1;
                                } else {
                                    pmap[npok - 1] = 0;
                                    npok--;
                                    incnpok = 0;
                                }
                            }
                        } else {
                            npok = rule->GetNOldP();
                            incnpok = 0;

                            if (ok)
                                foundmap[ri - 1]++;

                            ok = 1;

                            // check orientations

                            for (size_t i = 1; i <= rule->GetNOrientations(); i++) {
                                if (CW(lpoints.Get(pmap.Get(rule->GetOrientation(i).i1)),
                                       lpoints.Get(pmap.Get(rule->GetOrientation(i).i2)),
                                       lpoints.Get(pmap.Get(rule->GetOrientation(i).i3)))) {
                                    ok = 0;
                                    break;
                                }
                            }

                            if (!ok) continue;

                            Vector oldu(2 * rule->GetNOldP());

                            for (size_t i = 1; i <= rule->GetNOldP(); i++) {
                                Vec2d ui(rule->GetPoint(i), lpoints.Get(pmap.Get(i)));
                                oldu(2 * i - 2) = ui.X();
                                oldu(2 * i - 1) = ui.Y();
                            }

                            rule->SetFreeZoneTransformation(oldu, tolerance);

                            if (!rule->ConvexFreeZone()) {
                                ok = 0;
                            }

                            // check freezone:
                            if (!ok) continue;
                            for (int i = 1; i <= maxlegalpoint && ok; i++) {
                                if (!pused.Get(i) &&
                                    rule->IsInFreeZone(lpoints.Get(i))) {
                                    ok = 0;
                                    break;
                                }
                            }

                            if (!ok) continue;
                            for (int i = maxlegalpoint + 1; i <= lpoints.size(); i++) {
                                if (rule->IsInFreeZone(lpoints.Get(i))) {
                                    ok = 0;
                                    break;
                                }
                            }

                            if (!ok) continue;
                            for (int i = 1; i <= maxlegalline; i++) {
                                if (!lused.Get(i) &&
                                    rule->IsLineInFreeZone(lpoints.Get(llines.Get(i).I1()),
                                                           lpoints.Get(llines.Get(i).I2()))) {
                                    ok = 0;
                                    break;
                                }
                            }

                            if (!ok) continue;

                            for (int i = maxlegalline + 1; i <= llines.size(); i++) {
                                if (rule->IsLineInFreeZone(lpoints.Get(llines.Get(i).I1()),
                                                           lpoints.Get(llines.Get(i).I2()))) {
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
                                for (size_t i = oldnp + 1; i <= rule->GetNP(); i++) {
                                    Point2d np = rule->GetPoint(i);
                                    np.X() += newu(2 * (i - oldnp) - 2);
                                    np.Y() += newu(2 * (i - oldnp) - 1);

                                    pmap[i - 1] = lpoints.push_back(np);
                                }
                            }

                            // Setze neue Linien:
                            for (size_t i = rule->GetNOldL() + 1; i <= rule->GetNL(); i++) {
                                llines.push_back(INDEX_2(
                                        pmap.Get(rule->GetLine(i)[0]),
                                        pmap.Get(rule->GetLine(i)[1])));
                            }

                            // delete old lines:
                            for (size_t i = 1; i <= rule->GetNDelL(); i++) {
                                dellines.push_back(sortlines[lmap.Get(rule->GetDelLine(i)) - 1]);
                            }

                            // insert new elements:
                            for (size_t i = 0; i < rule->GetNE(); i++) {
                                elements.push_back(rule->GetElement(i + 1));
                                for (size_t j = 1; j <= elements[i].GetNP(); j++) {
                                    elements[i].PNum(j) = pmap.Get(elements[i].PNum(j));
                                }
                            }

                            double elerr = 0;
                            for (int i = 1; i <= elements.size(); i++) {
                                double hf;
                                if (!mp.quad) {
                                    hf = CalcElementBadness(lpoints, elements[i - 1]);
                                } else {
                                    // FIXME: this seems to be bugged
//                                    hf = elements[i - 1].CalcJacobianBadness(lpoints) * 5;
                                    hf = 1.0;
                                }

                                if (hf > elerr) elerr = hf;
                            }

                            canuse[ri - 1]++;

                            if (elerr < 0.99 * minelerr) {
                                minelerr = elerr;
                                found = ri;

                                tempnewpoints = lpoints.Range(noldlp, lpoints.size());
                                tempnewlines = llines.Range(noldll, llines.size());
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

                    nlok = rule->GetNOldL();

                    lused[lmap.Get(nlok)-1] = 0;

                    for (int j = 1; j <= 2; j++) {
                        int refpi = rule->GetPointNr(nlok, j);
                        pused[pmap.Get(refpi) - 1]--;

                        if (pused.Get(pmap.Get(refpi)) == 0)
                            pmap[refpi-1] = 0;
                    }
                }
            }
        }

        if (found) {
            lpoints.Append(tempnewpoints);
            llines1.Append(tempnewlines);
            dellines.Append(tempdellines);
            elements.Append(tempelements);
        }

        return found;
    }
}  // namespace meshit
