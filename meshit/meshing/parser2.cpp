#include <sstream>

#include "../meshit.hpp"
#include "meshing2.hpp"

#ifdef WIN32
#define COMMASIGN ':'
#else
#define COMMASIGN ','
#endif

namespace meshit {

    void LoadMatrixLine(std::istream & ist, DenseMatrix & m, int line)
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

            if (ch == 'x' || ch == 'X')
                m.Elem(line, 2 * pnum - 1) = f;
            if (ch == 'y' || ch == 'Y')
                m.Elem(line, 2 * pnum) = f;

            ist >> ch;
            if (ch == COMMASIGN)
                ist >> ch;
        }
    }

    void netrule::LoadRule(std::istream & ist)
    {
        char buf[256];
        char ch;
        Point2d p;
        INDEX_2 lin;
        int i, j;
        DenseMatrix tempoldutonewu(20, 20), tempoldutofreearea(20, 20),
                tempoldutofreearealimit(20, 20);

        tempoldutonewu = 0;
        tempoldutofreearea = 0;
        tempoldutofreearealimit = 0;

        noldp = 0;
        noldl = 0;

        ist.get(buf, sizeof (buf), '"');
        ist.get(ch);
        ist.get(buf, sizeof (buf), '"');
        ist.get(ch);

        // if(name != NULL) 
        delete [] name;
        name = new char[strlen(buf) + 1];
        strcpy(name, buf);

        do {
            ist >> buf;

            if (strcmp(buf, "quality") == 0) {
                ist >> quality;
            }

            else if (strcmp(buf, "mappoints") == 0) {
                ist >> ch;

                while (ch == '(') {
                    ist >> p.X();
                    ist >> ch; // ','
                    ist >> p.Y();
                    ist >> ch; // ')'

                    points.push_back(p);
                    noldp++;

                    tolerances.resize(noldp);
                    tolerances[noldp - 1].f1 = 1.0;
                    tolerances[noldp - 1].f2 = 0;
                    tolerances[noldp - 1].f3 = 1.0;

                    ist >> ch;
                    while (ch != ';') {
                        if (ch == '{') {
                            ist >> tolerances[noldp - 1].f1;
                            ist >> ch; // ','
                            ist >> tolerances[noldp - 1].f2;
                            ist >> ch; // ','
                            ist >> tolerances[noldp - 1].f3;
                            ist >> ch; // '}'
                        }
                        else if (ch == 'd') {
                            //            delpoints.Append (noldp);
                            ist >> ch; // 'e'
                            ist >> ch; // 'l'
                        }

                        ist >> ch;
                    }

                    ist >> ch;
                }

                ist.putback(ch);
            }


            else if (strcmp(buf, "maplines") == 0) {
                ist >> ch;

                while (ch == '(') {
                    ist >> lin.I1();
                    ist >> ch; // ','
                    ist >> lin.I2();
                    ist >> ch; // ')'

                    lines.push_back(lin);
                    linevecs.push_back(points[lin.I2() - 1] - points[lin.I1() - 1]);
                    noldl++;
                    linetolerances.resize(noldl);
                    linetolerances[noldl - 1].f1 = 0;
                    linetolerances[noldl - 1].f2 = 0;
                    linetolerances[noldl - 1].f3 = 0;

                    ist >> ch;
                    while (ch != ';') {
                        if (ch == '{') {
                            ist >> linetolerances[noldl - 1].f1;
                            ist >> ch; // ','
                            ist >> linetolerances[noldl - 1].f2;
                            ist >> ch; // ','
                            ist >> linetolerances[noldl - 1].f3;
                            ist >> ch; // '}'
                        }
                        else if (ch == 'd') {
                            dellines.push_back(noldl);
                            ist >> ch; // 'e'
                            ist >> ch; // 'l'
                        }

                        ist >> ch;
                    }

                    ist >> ch;
                }

                ist.putback(ch);
            }

            else if (strcmp(buf, "newpoints") == 0) {
                ist >> ch;

                while (ch == '(') {
                    ist >> p.X();
                    ist >> ch; // ','
                    ist >> p.Y();
                    ist >> ch; // ')'

                    points.push_back(p);

                    ist >> ch;
                    while (ch != ';') {
                        if (ch == '{') {
                            LoadMatrixLine(ist, tempoldutonewu,
                                    2 * (points.size() - noldp) - 1);

                            ist >> ch; // '{'
                            LoadMatrixLine(ist, tempoldutonewu,
                                    2 * (points.size() - noldp));
                        }

                        ist >> ch;
                    }

                    ist >> ch;
                }

                ist.putback(ch);
            }

            else if (strcmp(buf, "newlines") == 0) {
                ist >> ch;

                while (ch == '(') {
                    ist >> lin.I1();
                    ist >> ch; // ','
                    ist >> lin.I2();
                    ist >> ch; // ')'

                    lines.push_back(lin);
                    linevecs.push_back(points[lin.I2() - 1] - points[lin.I1() - 1]);

                    ist >> ch;
                    while (ch != ';') {
                        ist >> ch;
                    }

                    ist >> ch;
                }

                ist.putback(ch);
            }

            else if (strcmp(buf, "freearea") == 0) {
                ist >> ch;

                while (ch == '(') {
                    ist >> p.X();
                    ist >> ch; // ','
                    ist >> p.Y();
                    ist >> ch; // ')'

                    freezone.push_back(p);
                    freezonelimit.push_back(p);

                    ist >> ch;
                    while (ch != ';') {
                        if (ch == '{') {
                            LoadMatrixLine(ist, tempoldutofreearea,
                                    2 * freezone.size() - 1);

                            ist >> ch; // '{'
                            LoadMatrixLine(ist, tempoldutofreearea,
                                    2 * freezone.size());
                        }

                        ist >> ch;
                    }

                    ist >> ch;
                }

                for (i = 1; i <= tempoldutofreearealimit.Height(); i++)
                    for (j = 1; j <= tempoldutofreearealimit.Width(); j++)
                        tempoldutofreearealimit.Elem(i, j) =
                            tempoldutofreearea.Elem(i, j);


                ist.putback(ch);
            }
            else if (strcmp(buf, "freearea2") == 0) {
                ist >> ch;
                int freepi = 0;
                tempoldutofreearealimit = 0;

                while (ch == '(') {
                    freepi++;

                    ist >> p.X();
                    ist >> ch; // ','
                    ist >> p.Y();
                    ist >> ch; // ')'

                    freezonelimit[freepi - 1] = p;

                    ist >> ch;
                    while (ch != ';') {
                        if (ch == '{') {
                            LoadMatrixLine(ist, tempoldutofreearealimit,
                                    2 * freepi - 1);

                            ist >> ch; // '{'
                            LoadMatrixLine(ist, tempoldutofreearealimit,
                                    2 * freepi);
                        }

                        ist >> ch;
                    }

                    ist >> ch;
                }

                ist.putback(ch);
            }

            else if (strcmp(buf, "elements") == 0) {
                ist >> ch;

                while (ch == '(') {
                    elements.push_back(Element2d());
                    //                    Element2d & last = *elements.end();
                    Element2d & last = elements[elements.size() - 1];

                    ist >> last.PNum(1);
                    ist >> ch; // ','

                    if (ch == COMMASIGN) {
                        ist >> last.PNum(2);
                        ist >> ch; // ','
                    }
                    if (ch == COMMASIGN) {
                        ist >> last.PNum(3);
                        ist >> ch; // ','
                    }
                    if (ch == COMMASIGN) {
                        last.SetType(QUAD);
                        ist >> last.PNum(4);
                        ist >> ch; // ','
                    }

                    ist >> ch;
                    while (ch != ';') {
                        ist >> ch;
                    }

                    ist >> ch;
                }

                ist.putback(ch);
            }

            else if (strcmp(buf, "orientations") == 0) {
                ist >> ch;

                while (ch == '(') {
                    //        threeint a = threeint();
                    orientations.push_back(threeint());
                    //                    threeint & last = *orientations.end();
                    threeint & last = orientations[orientations.size() - 1];

                    ist >> last.i1;
                    ist >> ch; // ','
                    ist >> last.i2;
                    ist >> ch; // ','
                    ist >> last.i3;
                    ist >> ch; // ','

                    ist >> ch;
                    while (ch != ';') {
                        ist >> ch;
                    }

                    ist >> ch;
                }

                ist.putback(ch);
            }

            else if (strcmp(buf, "endrule") != 0) {
                MESHIT_LOG_ERROR("Parser error, unknown token " << buf);
            }
        } while (!ist.eof() && strcmp(buf, "endrule") != 0);

        oldutonewu.SetSize(2 * (points.size() - noldp), 2 * noldp);
        oldutofreearea.SetSize(2 * freezone.size(), 2 * noldp);
        oldutofreearealimit.SetSize(2 * freezone.size(), 2 * noldp);

        for (i = 1; i <= oldutonewu.Height(); i++)
            for (j = 1; j <= oldutonewu.Width(); j++)
                oldutonewu.Elem(i, j) = tempoldutonewu.Elem(i, j);

        for (i = 1; i <= oldutofreearea.Height(); i++)
            for (j = 1; j <= oldutofreearea.Width(); j++)
                oldutofreearea.Elem(i, j) = tempoldutofreearea.Elem(i, j);

        for (i = 1; i <= oldutofreearea.Height(); i++)
            for (j = 1; j <= oldutofreearea.Width(); j++)
                oldutofreearealimit.Elem(i, j) = tempoldutofreearealimit.Elem(i, j);

        freesetinequ.SetSize(freezone.size());


        {
            char ok;
            int minn;
            Array<int> pnearness(noldp);

            for (i = 1; i <= pnearness.size(); i++)
                pnearness.Elem(i) = 1000;

            for (j = 1; j <= 2; j++)
                pnearness.Elem(GetPointNr(1, j)) = 0;

            do {
                ok = 1;

                for (i = 1; i <= (int) noldl; i++) {
                    minn = 1000;
                    for (j = 1; j <= 2; j++)
                        minn = std::min(minn, pnearness.Get(GetPointNr(i, j)));

                    for (j = 1; j <= 2; j++)
                        if (pnearness.Get(GetPointNr(i, j)) > minn + 1) {
                            ok = 0;
                            pnearness.Elem(GetPointNr(i, j)) = minn + 1;
                        }
                }
            } while (!ok);

            lnearness.resize(noldl);

            for (size_t i = 0; i < noldl; i++) {
                lnearness[i] = 0;
                for (j = 1; j <= 2; j++) {
                    lnearness[i] += pnearness.Get(GetPointNr(i + 1, j));
                }
            }
        }

        oldutofreearea_i.resize(10);
        freezone_i.resize(10);

        for (size_t i = 0; i < oldutofreearea_i.size(); i++) {
            double lam1 = 1.0 / (i + 1);

            oldutofreearea_i[i] = new DenseMatrix(oldutofreearea.Height(), oldutofreearea.Width());
            DenseMatrix & mati = *oldutofreearea_i[i];
            for (j = 0; j < oldutofreearea.Height(); j++)
                for (size_t k = 0; k < oldutofreearea.Width(); k++)
                    mati(j, k) = lam1 * oldutofreearea(j, k) + (1 - lam1) * oldutofreearealimit(j, k);

            freezone_i[i] = new Array<Point2d> (freezone.size());
            Array<Point2d> & fzi = *freezone_i[i];
            for (size_t j = 0; j < freezone.size(); j++)
                fzi[j] = freezonelimit[j] + lam1 * (freezone[j] - freezonelimit[j]);
        }
    }

    extern const char * triarules[];
    extern const char * quadrules[];

    void Meshing2::LoadRules(const char * filename, bool quad)
    {
        char buf[256];
        std::istream * ist;
        std::string tr1;

        if (filename) {
            ist = new std::ifstream(filename);
        }
        else {
            /* connect tetrules to one string */
            const char ** hcp;

            if (!quad) {
                hcp = triarules;
                MESHIT_LOG_DEBUG("load internal triangle rules");
            }
            else {
                hcp = quadrules;
                MESHIT_LOG_DEBUG("load internal quad rules");
            }

            size_t len = 0;
            while (*hcp) {
                len += strlen(*hcp);
                hcp++;
            }
            tr1.reserve(len + 1);

            if (!quad)
                hcp = triarules;
            else
                hcp = quadrules;

            while (*hcp) {
                tr1.append(*hcp);
                hcp++;
            }

#ifdef WIN32
            // VC++ 2005 workaround
            for (std::string::size_type i = 0; i < tr1.size(); i++) {
                if (tr1[i] == ',')
                    tr1[i] = ':';
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
                netrule * rule = new netrule;
                rule->LoadRule(*ist);
                rules.push_back(rule);
            }
        }
        delete ist;
    }

}
