#pragma once

namespace steps::model {

class Chan;
class ChanState;
class Complex;
class ComplexReac;
class ComplexSReac;
class Diff;
class Endocytosis;
class Exocytosis;
class GHKcurr;
class LinkSpec;
class Model;
class OhmicCurr;
class Raft;
class RaftDis;
class RaftEndocytosis;
class RaftGen;
class RaftSReac;
class Raftsys;
class Reac;
class Spec;
class SReac;
class Surfsys;
class VDepSReac;
class VesBind;
class Vesicle;
class VesSDiff;
class VesSReac;
class VesSurfsys;
class VesUnbind;
class Volsys;

enum Immobilization { IMMOBILIZING = 1, NO_EFFECT = 0, MOBILIZING = -1 };

}  // namespace steps::model
