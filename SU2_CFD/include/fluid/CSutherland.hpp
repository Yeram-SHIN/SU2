/*!
 * \file CSutherland.hpp
 * \brief Defines Sutherland's Law for laminar viscosity.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna, T. Economon
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "CConductivityModel.hpp"

/*!
 * \class CSutherland
 * \brief Defines Sutherland's Law for laminar conductivity.
 * \author S.Vitale, M.Pini
 */
class CCondSutherland final : public CConductivityModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CSutherland(su2double mu_cond_ref, su2double t_cond_ref, su2double s_cond) : mu_cond_ref_(mu_cond_ref), t_cond_ref_(t_cond_ref), s_cond_(s_cond) {}

  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double t, su2double rho) override {
    const su2double t_ref_inv = 1.0 / t_cond_ref_;
    const su2double t_nondim = t_ref_inv * t;
    mu_ = mu_cond_ref_ * t_nondim * sqrt(t_nondim) * ((t_cond_ref_ + s_cond_) / (t + s_cond_));

    /*--- Set Viscosity Derivatives. ---*/
    const su2double ts_inv = 1.0 / (t + s_cond_);
    dmudrho_t_ = 0.0;
    dmudt_rho_ = mu_cond_ref_ * (t_cond_ref_ + s_cond_) * ts_inv * sqrt(t_nondim) * (1.5 * t_ref_inv - t_nondim * ts_inv);
  }

 private:
  const su2double mu_cond_ref_{0.0};    /*!< \brief Internal Energy. */
  const su2double t_cond_ref_{0.0};     /*!< \brief DpDd_e. */
  const su2double s_cond_{0.0};         /*!< \brief DpDe_d. */
};
