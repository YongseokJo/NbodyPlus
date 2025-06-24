#ifndef AR_PERTURBER
#define AR_PERTURBER
#ifdef FEWBODY

#include "Common/list.h"
#include "Common/binary_tree.h"
#include "../particle.h"

//! Perturber class for AR integration
class Perturber {
public:

    Float soft_pert_min; ///> minimum soft perturbation

    Perturber(): soft_pert_min(Float(0.0)) {}

    //! clear function
    void clear() {
    }

    //! check parameters status
    bool checkParams() {
        return true;
    }

    //! calculate soft_pert_min for slowdown pert_out
    /*! \Delta F = G m_cm m_p (apo) / rp^3
        Pert_out = \Delta F /(G apo)
        @param[in] _bin: binary parameters
        @param[in] _G: gravitatioal constant
     */
    // template <class Tptcl>
    void calcSoftPertMin(const AR::BinaryTree<Particle>& _bin, const Float _G) {
        soft_pert_min = 0.0;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ostream & _fout, const int _width=20){
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
    } 

};
#endif
#endif