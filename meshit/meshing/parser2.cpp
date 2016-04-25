#include <sstream>

#include "mesh_generator.hpp"

#ifdef WIN32
#define COMMASIGN ':'
#else
#define COMMASIGN ','
#endif

namespace meshit {

void LoadMatrixLine(std::istream& ist, DenseMatrix& m, int line)
{
    char ch;
    int pnum;
    float f;

    ist >> ch;
    while (ch != '}') {
        ist.putback(ch);
        ist >> f;
        ist >> ch;
        ist >> pnum;

        if (ch == 'x' || ch == 'X') m.Elem(line - 1, 2 * pnum - 2) = f;
        if (ch == 'y' || ch == 'Y') m.Elem(line - 1, 2 * pnum - 1) = f;

        ist >> ch;
        if (ch == COMMASIGN) ist >> ch;
    }
}

void netrule::LoadRule(std::istream& ist)
{
    char buf[256];
    char ch;
    Point2d p;
    DenseMatrix tempoldutonewu(20, 20);
    DenseMatrix tempoldutofreearea(20, 20);
    DenseMatrix tempoldutofreearealimit(20, 20);

    tempoldutonewu = 0;
    tempoldutofreearea = 0;
    tempoldutofreearealimit = 0;

    noldp = 0;
    noldl = 0;

    ist.get(buf, sizeof(buf), '"');
    ist.get(ch);
    ist.get(buf, sizeof(buf), '"');
    ist.get(ch);

    delete[] name;
    name = new char[strlen(buf) + 1];
    strcpy(name, buf);

    do {
        ist >> buf;

        if (strcmp(buf, "quality") == 0) {
            ist >> quality;
        } else if (strcmp(buf, "mappoints") == 0) {
            ist >> ch;
            while (ch == '(') {
                ist >> p.X() >> ch;  // ','
                ist >> p.Y() >> ch;  // ')'

                points.push_back(p);
                tolerances.push_back({1.0, 0.0, 1.0});
                noldp++;

                ist >> ch;
                while (ch != ';') {
                    if (ch == '{') {
                        ist >> tolerances.back().f1 >> ch;  // ','
                        ist >> tolerances.back().f2 >> ch;  // ','
                        ist >> tolerances.back().f3 >> ch;  // '}'
                    } else if (ch == 'd') {
                        ist >> ch;  // 'e'
                        ist >> ch;  // 'l'
                    }
                    ist >> ch;
                }
                ist >> ch;
            }
            ist.putback(ch);
        } else if (strcmp(buf, "maplines") == 0) {
            ist >> ch;
            while (ch == '(') {
                INDEX li1, li2;
                ist >> li1 >> ch;  // ','
                ist >> li2 >> ch;  // ')'
                li1--;  // indices are stored 0-based
                li2--;
                lines.push_back(INDEX_2(li1, li2));
                linevecs.push_back(points[li2] - points[li1]);
                linetolerances.push_back({0.0, 0.0, 0.0});
                noldl++;

                ist >> ch;
                while (ch != ';') {
                    if (ch == '{') {
                        ist >> linetolerances.back().f1 >> ch;  // ','
                        ist >> linetolerances.back().f2 >> ch;  // ','
                        ist >> linetolerances.back().f3 >> ch;  // '}'
                    } else if (ch == 'd') {
                        dellines.push_back(noldl);
                        ist >> ch;  // 'e'
                        ist >> ch;  // 'l'
                    }
                    ist >> ch;
                }
                ist >> ch;
            }
            ist.putback(ch);
        } else if (strcmp(buf, "newpoints") == 0) {
            ist >> ch;
            while (ch == '(') {
                ist >> p.X() >> ch;  // ','
                ist >> p.Y() >> ch;  // ')'
                points.push_back(p);
                ist >> ch;
                while (ch != ';') {
                    if (ch == '{') {
                        LoadMatrixLine(ist, tempoldutonewu, 2 * (points.size() - noldp) - 1);
                        ist >> ch;  // '{'
                        LoadMatrixLine(ist, tempoldutonewu, 2 * (points.size() - noldp));
                    }
                    ist >> ch;
                }
                ist >> ch;
            }
            ist.putback(ch);
        } else if (strcmp(buf, "newlines") == 0) {
            ist >> ch;
            while (ch == '(') {
                INDEX li1, li2;
                ist >> li1 >> ch;  // ','
                ist >> li2 >> ch;  // ')'
                li1--;  // indices are stored 0-based
                li2--;
                lines.push_back(INDEX_2(li1, li2));
                linevecs.push_back(points[li2] - points[li1]);
                ist >> ch;
                while (ch != ';') {
                    ist >> ch;
                }
                ist >> ch;
            }
            ist.putback(ch);
        } else if (strcmp(buf, "freearea") == 0) {
            ist >> ch;
            while (ch == '(') {
                ist >> p.X() >> ch;  // ','
                ist >> p.Y() >> ch;  // ')'
                freezone.push_back(p);
                freezonelimit.push_back(p);
                ist >> ch;
                while (ch != ';') {
                    if (ch == '{') {
                        LoadMatrixLine(ist, tempoldutofreearea, 2 * freezone.size() - 1);
                        ist >> ch;  // '{'
                        LoadMatrixLine(ist, tempoldutofreearea, 2 * freezone.size());
                    }
                    ist >> ch;
                }
                ist >> ch;
            }
            for (size_t i = 0; i < tempoldutofreearealimit.Height(); i++) {
                for (size_t j = 0; j < tempoldutofreearealimit.Width(); j++) {
                    tempoldutofreearealimit.Elem(i, j) = tempoldutofreearea.Elem(i, j);
                }
            }
            ist.putback(ch);
        } else if (strcmp(buf, "freearea2") == 0) {
            ist >> ch;
            int freepi = 0;
            tempoldutofreearealimit = 0;
            while (ch == '(') {
                freepi++;
                ist >> p.X() >> ch;  // ','
                ist >> p.Y() >> ch;  // ')'
                freezonelimit[freepi - 1] = p;
                ist >> ch;
                while (ch != ';') {
                    if (ch == '{') {
                        LoadMatrixLine(ist, tempoldutofreearealimit, 2 * freepi - 1);
                        ist >> ch;  // '{'
                        LoadMatrixLine(ist, tempoldutofreearealimit, 2 * freepi);
                    }
                    ist >> ch;
                }
                ist >> ch;
            }
            ist.putback(ch);
        } else if (strcmp(buf, "elements") == 0) {
            ist >> ch;
            while (ch == '(') {
                elements.push_back(Element2d());
                Element2d& last = elements.back();

                ist >> last.PointID(0);
                ist >> ch;  // ','

                if (ch == COMMASIGN) {
                    ist >> last.PointID(1) >> ch;  // ','
                }
                if (ch == COMMASIGN) {
                    ist >> last.PointID(2) >> ch;  // ','
                }
                if (ch == COMMASIGN) {
                    MESHIT_LOG_FATAL("netrule::LoadRule has just read a QUAD elements. Not supported anymore.");
                    exit(1);
                }
                ist >> ch;
                while (ch != ';') {
                    ist >> ch;
                }
                ist >> ch;
            }
            ist.putback(ch);
        } else if (strcmp(buf, "orientations") == 0) {
            ist >> ch;
            while (ch == '(') {
                int i1, i2, i3;
                ist >> i1 >> ch;  // ','
                ist >> i2 >> ch;  // ','
                ist >> i3 >> ch;  // ','
                ist >> ch;
                orientations.push_back(threeint{i1, i2, i3});
                while (ch != ';') {
                    ist >> ch;
                }
                ist >> ch;
            }
            ist.putback(ch);
        } else if (strcmp(buf, "endrule") != 0) {
            MESHIT_LOG_ERROR("Parser error, unknown token " << buf);
        }
    } while (!ist.eof() && strcmp(buf, "endrule") != 0);

    oldutonewu.SetSize(2 * (points.size() - noldp), 2 * noldp);
    oldutofreearea.SetSize(2 * freezone.size(), 2 * noldp);
    oldutofreearealimit.SetSize(2 * freezone.size(), 2 * noldp);

    for (size_t i = 0; i < oldutonewu.Height(); i++) {
        for (size_t j = 0; j < oldutonewu.Width(); j++) {
            oldutonewu.Elem(i, j) = tempoldutonewu.Elem(i, j);
        }
    }
    for (size_t i = 0; i < oldutofreearea.Height(); i++) {
        for (size_t j = 0; j < oldutofreearea.Width(); j++) {
            oldutofreearea.Elem(i, j) = tempoldutofreearea.Elem(i, j);
        }
    }
    for (size_t i = 0; i < oldutofreearea.Height(); i++) {
        for (size_t j = 0; j < oldutofreearea.Width(); j++) {
            oldutofreearealimit.Elem(i, j) = tempoldutofreearealimit.Elem(i, j);
        }
    }

    freesetinequ.SetSize(freezone.size());

    std::vector<uint32_t> pnearness(noldp, 1000);

    pnearness[lines[0].I1()] = 0;
    pnearness[lines[0].I2()] = 0;

    bool ok;
    do {
        ok = true;
        for (size_t i = 0; i < noldl; i++) {
            INDEX idx1 = lines[i].I1();
            INDEX idx2 = lines[i].I2();
            uint32_t minn = 1000;
            minn = std::min(minn, pnearness[idx1]);
            minn = std::min(minn, pnearness[idx2]);
            if (pnearness[idx1] > minn + 1) {
                pnearness[idx1] = minn + 1;
                ok = false;
            }
            if (pnearness[idx2] > minn + 1) {
                pnearness[idx2] = minn + 1;
                ok = false;
            }
        }
    } while (!ok);

    lnearness.resize(noldl);
    for (size_t i = 0; i < noldl; i++) {
        lnearness[i] = pnearness[lines[i].I1()] + pnearness[lines[i].I2()];
    }

    freezone_i.resize(10);
    oldutofreearea_i.resize(10);

    for (size_t i = 0; i < oldutofreearea_i.size(); i++) {
        double lam1 = 1.0 / (i + 1);

        oldutofreearea_i[i] = new DenseMatrix(oldutofreearea.Height(), oldutofreearea.Width());
        DenseMatrix& mati = *oldutofreearea_i[i];
        for (size_t j = 0; j < oldutofreearea.Height(); j++)
            for (size_t k = 0; k < oldutofreearea.Width(); k++)
                mati(j, k) = lam1 * oldutofreearea(j, k) + (1 - lam1) * oldutofreearealimit(j, k);

        freezone_i[i] = new std::vector<Point2d>(freezone.size());
        std::vector<Point2d>& fzi = *freezone_i[i];
        for (size_t j = 0; j < freezone.size(); j++) {
            fzi[j] = freezonelimit[j] + lam1 * (freezone[j] - freezonelimit[j]);
        }
    }
}

extern const char* triarules[];

void MeshGenerator::LoadRules(const char* filename)
{
    char buf[256];
    std::istream* ist;
    std::string tr1;

    if (filename) {
        ist = new std::ifstream(filename);
    } else {
        const char** hcp = triarules;
        MESHIT_LOG_DEBUG("load internal triangle rules");

        size_t len = 0;
        while (*hcp) {
            len += strlen(*hcp);
            hcp++;
        }
        tr1.reserve(len + 1);

        hcp = triarules;
        while (*hcp) {
            tr1.append(*hcp);
            hcp++;
        }

#ifdef WIN32
        // VC++ 2005 workaround
        for (std::string::size_type i = 0; i < tr1.size(); i++) {
            if (tr1[i] == ',') tr1[i] = ':';
        }
#endif
        ist = new std::istringstream(tr1);
    }

    if (!ist->good()) {
        std::cerr << "Rule description file " << filename << " not found" << std::endl;
        delete ist;
        exit(1);
    }

    while (!ist->eof()) {
        buf[0] = 0;
        (*ist) >> buf;

        if (strcmp(buf, "rule") == 0) {
            netrule* rule = new netrule;
            rule->LoadRule(*ist);
            rules.push_back(rule);
        }
    }
    delete ist;
}

}  // namespace meshit
