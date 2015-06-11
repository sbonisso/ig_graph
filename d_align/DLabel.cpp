#include "DLabel.hpp"

DLabel::DLabel() {
    label = "";
    score = -1;
    ident = -1;
}

DLabel::DLabel(std::string s, int sc, int iden) {
    label = s;
    score = sc;
    ident = iden;
}

DLabel::~DLabel() {}

bool DLabel::operator > (const DLabel& dl) const {
    return (score > dl.score);
}

bool DLabel::operator < (const DLabel& dl) const {
    return (score < dl.score);
}

ostream& operator<<(ostream& out, const DLabel& dl) {
    out<<"["<<dl.label<<", "<<dl.score<<", "<<dl.ident<<"]";
    return out;
}
