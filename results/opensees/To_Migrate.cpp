//https://github.com/peer-open-source/xara/blob/stable/SRC/element/Shell/ASDShellT3.cpp
//line 307
computeBMatrix

        // geometric data
        const auto& p1 = LCS.P1();
        const auto& p2 = LCS.P2();
        const auto& p3 = LCS.P3();
        double y12 = p1[1] - p2[1];
        double y23 = p2[1] - p3[1];
        double y31 = p3[1] - p1[1];
        double x23 = p2[0] - p3[0];
        double x31 = p3[0] - p1[0];
        double x12 = p1[0] - p2[0];
        double x32 = p3[0] - p2[0];
        double y32 = p3[1] - p2[1];
        double x13 = p1[0] - p3[0];
        double y13 = p1[1] - p3[1];
        double x21 = p2[0] - p1[0];
        double y21 = p2[1] - p1[1];

        // membrane part (ANDeS, basic)
        constexpr double alpha_membrane = 1.5;
        double A2 = LCS.Area() * 2.0;
        B(0, 0) = y23 / A2;
        B(2, 0) = -x23 / A2;
        B(1, 1) = -x23 / A2;
        B(2, 1) = y23 / A2;
        B(0, 5) = alpha_membrane * y23 * (-y31 + y12) / 6 / A2;
        B(1, 5) = alpha_membrane * (-x23) * (x31 - x12) / 6 / A2;
        B(2, 5) = alpha_membrane * (-x31 * y31 + x12 * y12) / 3 / A2;

        B(0, 6) = y31 / A2;
        B(2, 6) = -x31 / A2;
        B(1, 7) = -x31 / A2;
        B(2, 7) = y31 / A2;
        B(0, 11) = alpha_membrane * y31 * (-y12 + y23) / 6 / A2;
        B(1, 11) = alpha_membrane * (-x31) * (x12 - x23) / 6 / A2;
        B(2, 11) = alpha_membrane * (-x12 * y12 + x23 * y23) / 3 / A2;

        B(0, 12) = y12 / A2;
        B(2, 12) = -x12 / A2;
        B(1, 13) = -x12 / A2;
        B(2, 13) = y12 / A2;
        B(0, 17) = alpha_membrane * y12 * (-y23 + y31) / 6 / A2;
        B(1, 17) = alpha_membrane * (-x12) * (x23 - x31) / 6 / A2;
        B(2, 17) = alpha_membrane * (-x23 * y23 + x31 * y31) / 3 / A2;


        // bending part
        for (int i = 0; i < 3; ++i) {
            int ii = i * 6;
            B(3, ii + 4) = -dNdX(i, 0);
            B(4, ii + 3) = dNdX(i, 1);
            B(5, ii + 3) = dNdX(i, 0);
            B(5, ii + 4) = -dNdX(i, 1);
        }

        // shear part (MITC3 treatment of transverse shear locking)
        double phi1 = std::atan2(y21, x21);
        double phi2 = 3.141592653589793 * 0.5 - std::atan2(x31, y31);
        auto& BsO = ASDShellT3Globals::instance().BsO;
        BsO(0, 0) =  std::sin(phi2);
        BsO(0, 1) = -std::sin(phi1);
        BsO(1, 0) = -std::cos(phi2);
        BsO(1, 1) =  std::cos(phi1);
        auto& BsC = ASDShellT3Globals::instance().BsC;
        BsC(0, 0) = std::sqrt(x13 * x13 + y13 * y13) / A2;
        BsC(0, 1) = 0.0;
        BsC(1, 0) = 0.0;
        BsC(1, 1) = std::sqrt(x21 * x21 + y21 * y21) / A2;
        auto& BsT = ASDShellT3Globals::instance().BsT;
        BsT(0, 0) = -1.0;
        BsT(0, 1) = -(y21 + y32 * xi) / 2.0;
        BsT(0, 2) = (x21 + x32 * xi) / 2.0;
        BsT(0, 3) = 1.0;
        BsT(0, 4) = -(y21 + y13 * xi) / 2.0;
        BsT(0, 5) = (x21 + x13 * xi) / 2.0;
        BsT(0, 7) = -(y21 * xi) / 2.0;
        BsT(0, 8) = (x21 * xi) / 2.0;
        BsT(1, 0) = -1.0;
        BsT(1, 1) = (y13 + y32 * eta) / 2.0;
        BsT(1, 2) = -(x13 + x32 * eta) / 2.0;
        BsT(1, 4) = (y13 * eta) / 2.0;
        BsT(1, 5) = -(x13 * eta) / 2.0;
        BsT(1, 6) = 1.0;
        BsT(1, 7) = (y13 + y21 * eta) / 2.0;
        BsT(1, 8) = -(x13 + x21 * eta) / 2.0;
        // Bs (modified shear B matrix = O*C*T)
        auto& BsOC = ASDShellT3Globals::instance().BsOC;
        BsOC.addMatrixProduct(0.0, BsO, BsC, 1.0);
        auto& Bs = ASDShellT3Globals::instance().Bs;
        Bs.addMatrixProduct(0.0, BsOC, BsT, 1.0);
        // add it to B
        for (int i = 0; i < 2; ++i) {
            for (int node = 0; node < 3; ++node) {
                int j = node * 6;
                int k = node * 3;
                for (int dof = 0; dof < 3; ++dof) {
                    B(i + 6, j + dof + 2) = Bs(i, k + dof);
                }
            }
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    for (int igauss = 0; igauss < gx.size(); ++igauss)
    {
        // Current integration point data
        double xi = gx[igauss];
        double eta = gy[igauss];
        double w = gw[igauss];
        shapeFunctions(xi, eta, N);
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);
        double dA = w * jac.detJ;
        dNdX.addMatrixProduct(0.0, dN, jac.invJ, 1.0);

        // Strain-displacement matrix
        computeBMatrix(reference_cs, dNdX, N, xi, eta, !m_reduced_integration, B);
        if (m_reduced_integration) {
            computeBdrilling(dNdX, N, Bd, Bhx, Bhy);
        }

        // Update strain strain
        if (options & OPT_UPDATE)
        {
            // Section deformation
            if (m_angle != 0.0) {
                auto& Elocal = ASDShellT3Globals::instance().Elocal;
                Elocal.addMatrixVector(0.0, B, UL, 1.0);
                E.addMatrixVector(0.0, Re, Elocal, 1.0);
            }
            else {
                E.addMatrixVector(0.0, B, UL, 1.0);
            }

            // apply Stenberg stabilization
            E(6) *= sh_stab;
            E(7) *= sh_stab;

            // Update section
            result += m_sections[igauss]->setTrialSectionDeformation(E);

            // Drilling strain Ed = Bd*UL
            if (m_reduced_integration) {
                m_drill_strain = 0.0;
                m_drill_curvature_x = 0.0;
                m_drill_curvature_y = 0.0;
                for (int i = 0; i < 18; i++) {
                    m_drill_strain += Bd(i) * UL(i);
                    m_drill_curvature_x += Bhx(i) * UL(i);
                    m_drill_curvature_y += Bhy(i) * UL(i);
                }
            }
        }

        // Invert bending terms for correct statement of equilibrium
        if ((options & OPT_RHS) || (options & OPT_LHS))
            invertBBendingTerms(B, B1);

        // Integrate RHS
        if (options & OPT_RHS)
        {
            // Section force
            if (m_angle != 0.0) {
                auto& Ssection = m_sections[igauss]->getStressResultant();
                S.addMatrixVector(0.0, Rs, Ssection, 1.0);
#ifdef OPS_USE_DAMPING
                if (m_damping[igauss]) {
                    m_damping[igauss]->update(Ssection);
                    auto& Sdsection = m_damping[igauss]->getDampingForce();
                    S.addMatrixVector(1.0, Rs, Sdsection, 1.0);
                }
#endif
            }
            else {
#ifdef OPS_USE_DAMPING
                S = m_sections[igauss]->getStressResultant();
                if (m_damping[igauss]) {
                    m_damping[igauss]->update(S);
                    S += m_damping[igauss]->getDampingForce();
                }
#endif
            }

            // apply Stenberg stabilization
            S(6) *= sh_stab;
            S(7) *= sh_stab;

            // Add current integration point contribution (RHS)
            RHS.addMatrixTransposeVector(1.0, B1, S, dA);

            // Drilling
            if (m_reduced_integration) {
                // Compute drilling damages
                if (m_drill_mode == DrillingDOF_NonLinear) {
                    drill_dstrain = m_sections[igauss]->getSectionDeformation();
                    drill_dstrain.addVector(1.0, m_nldrill->strain_comm, -1.0);
                    if (drill_dstrain.Norm() > 1.0e-10) {
                        drill_dstress = m_sections[igauss]->getStressResultant();
                        drill_dstress.addVector(1.0, m_nldrill->stress_comm, -1.0);
                        const Matrix& C0 = m_sections[igauss]->getInitialTangent();
                        drill_dstress_el.addMatrixVector(0.0, C0, drill_dstrain, 1.0);
                        double drill_damage = m_nldrill->damage_comm;
                        for (int j = 0; j < 8; ++j) {
                            double dx = drill_dstrain(j);
                            if (dx > 0.0) {
                                double dy = drill_dstress(j);
                                double dy0 = drill_dstress_el(j);
                                drill_damage = std::max(drill_damage, 1.0 - dy / dy0);
                            }
                        }
                        drill_damage = std::max(0.0, std::min(0.999, drill_damage));
                        m_nldrill->damage = drill_damage;
                    }
                }

                // Add drilling stress = Bd'*Sd * dA (RHS)
                double Sd = m_drill_stiffness * m_drill_strain;
                double Shx = m_drill_stiffness * m_drill_curvature_x * dh_scale_factor;
                double Shy = m_drill_stiffness * m_drill_curvature_y * dh_scale_factor;
                if (m_drill_mode == DrillingDOF_NonLinear) {
                    Sd *= (1.0 - m_nldrill->damage_comm);
                    Shx *= (1.0 - m_nldrill->damage_comm);
                    Shy *= (1.0 - m_nldrill->damage_comm);
                }
                for (int i = 0; i < 18; i++) {
                    RHS(i) += Bd(i) * Sd * dA;
                    RHS(i) += Bhx(i) * Shx * dA;
                    RHS(i) += Bhy(i) * Shy * dA;
                }
            }
        }

        // Integrate LHS
        if (options & OPT_LHS)
        {
            // Section tangent
            if (m_angle != 0.0) {
                Dsection = (options & OPT_LHS_IS_INITIAL) ?
                    m_sections[igauss]->getInitialTangent() :
                    m_sections[igauss]->getSectionTangent();
#ifdef OPS_USE_DAMPING
                if (m_damping[igauss]) 
                    Dsection *= m_damping[igauss]->getStiffnessMultiplier();
#endif
                auto& RsT = ASDShellT3Globals::instance().RsT;
                RsT.addMatrixTranspose(0.0, Rs, 1.0);
                auto& DRsT = ASDShellT3Globals::instance().DRsT;
                DRsT.addMatrixProduct(0.0, Dsection, RsT, 1.0);
                D.addMatrixProduct(0.0, Rs, DRsT, 1.0);
            }
            else {
                D = (options & OPT_LHS_IS_INITIAL) ?
                    m_sections[igauss]->getInitialTangent() :
                    m_sections[igauss]->getSectionTangent();
#ifdef OPS_USE_DAMPING
                if (m_damping[igauss])
                    D *= m_damping[igauss]->getStiffnessMultiplier();
#endif
            }

            // apply Stenberg stabilization
            for (int i = 0; i < 8; ++i) {
                D(6, i) *= sh_stab;
                D(7, i) *= sh_stab;
                D(i, 6) *= sh_stab;
                D(i, 7) *= sh_stab;
            }

            // Add current integration point contribution (LHS)
            B1TD.addMatrixTransposeProduct(0.0, B1, D, dA);
            LHS.addMatrixProduct(1.0, B1TD, B, 1.0);

            // Add drilling stiffness = Bd'*Kd*Bd * dA (LHS)
            if (m_reduced_integration) {
                double drill_tang = m_drill_stiffness;
                if (m_drill_mode == DrillingDOF_NonLinear)
                    drill_tang *= (1.0 - m_nldrill->damage_comm);
                for (int i = 0; i < 18; i++) {
                    for (int j = 0; j < 18; j++) {
                        LHS(i, j) += Bd(i) * drill_tang * Bd(j) * dA;
                        LHS(i, j) += Bhx(i) * drill_tang * Bhx(j) * dA * dh_scale_factor;
                        LHS(i, j) += Bhy(i) * drill_tang * Bhy(j) * dA * dh_scale_factor;
                    }
                }
            }
        }
    }

    // Transform LHS to global coordinate system
    m_transformation->transformToGlobal(local_cs, UG, UL, LHS, RHS, (options & OPT_LHS));

    // Subtract external loads if any
    if ((options & OPT_RHS) && m_load)
        RHS.addVector(1.0, *m_load, -1.0);

    // Done
    return result;