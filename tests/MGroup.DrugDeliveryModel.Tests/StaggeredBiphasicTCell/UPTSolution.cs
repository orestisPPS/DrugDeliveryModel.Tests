using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using TriangleNet.Meshing.Algorithm;
using System.Xml.Linq;
using MGroup.DrugDeliveryModel.Tests.PreliminaryModels;
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;
using BC = MGroup.DrugDeliveryModel.Tests.Commons.BoundaryAndInitialConditionsUtility.BoundaryConditionCase;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class UPTSolution
    {
        const double Sc = 0.1;

        private const double timeStep = 1E-5; // in sec
        const double totalTime = 10E-5; // in sec
        //const double totalTime = 0.0001 ; // in sec
        static int incrementsPertimeStep = 1;
        static int currentTimeStep = 0;

        #region Structural model properties

        static double density = 1;
        static double miNormal = 5; //KPa
        static double kappaNormal = 6.667; //Kpa

        static double miTumor = 22.44; //Kpa
        static double kappaTumor = 216.7; //Kpa

        static double lambda0 = 1;
        int nGaussPoints = 1;
        static double initial_dp_dx = 0.0;
        static double initial_dp_dy = 0.0;
        static double initial_dp_dz = 0.0;
        static double velocityDivInitialVal = 0;
        static Dictionary<double, double[]> Solution = new Dictionary<double, double[]>();
        private static List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>();

        #endregion

        #region Structural model BCs and Loads

        private static List<(BC, StructuralDof[], double[][], double[])> structuralDirichletBC =
            new List<(BC, StructuralDof[], double[][], double[])>()
                {
                    (BC.BottomDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationZ }, new double[1][]{new double[3] {0,0,0}}, new double[] { 0.0 }),
                    (BC.LeftDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationX }, new double[1][]{new double[3] {0,0,0}}, new double[] { 0.0 }),
                    (BC.RightDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationX }, new double[1][]{new double[3] {0.1,0,0}}, new double[] { 0.0 }),
                    (BC.FrontDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationY }, new double[1][]{new double[3] {0,0,0}}, new double[] { 0.0 }),
                    (BC.BackDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationY }, new double[1][]{new double[3] {0,0.1,0}}, new double[] { 0.0 }),

                };

        static double[] coordsLoad1 = new double[3] { 0.0, 0.0, 0.1 };
        static double[] coordsLoad2 = new double[3] { 0.1, 0.0, 0.1 };
        static double[] coordsLoad3 = new double[3] { 0.0, 0.1, 0.1 };
        static double[] coordsLoad4 = new double[3] { 0.1, 0.1, 0.1 };
        static double[][] loadCoords = new double[4][] { coordsLoad1, coordsLoad2, coordsLoad3, coordsLoad4 };
        static List<(BC, StructuralDof[], double[][], double[])> structuralNeumannBC =
            new List<(BC, StructuralDof[], double[][], double[])>()
            {
                (BC.TopPointFlux, new StructuralDof[1]{StructuralDof.TranslationZ}, loadCoords, new double []{-1E-4 / 4d})
            };

        #endregion

        #region ToDo Orestis log task 1 
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: logged Doftype, 4: found node Id, 5: results 
        /// </summary>
        static List<(double[], string, StructuralDof, int, double[])> nodeDisplacementLogs = new List<(double[], string, StructuralDof, int, double[])>()
            {(new double[]{ 0.04930793848882013,0.04994681648346263,0.075 }, "CornerNodeTranslationZ.txt",StructuralDof.TranslationZ,-1, new double[0])};

        static double[] structuralMonitorNodeCoords = new double[] { 0.0, 0.0, 0.1 };
        private static int structuralMonitorID;
        static StructuralDof structuralMonitorDOF = StructuralDof.TranslationZ;

        #endregion

        #region ToDo Orestis log task 2
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: found element Id, 4: results 
        /// </summary>
        static List<(double[], string, int, double[][])> gpVelocityGraadientLogs = new List<(double[], string, int, double[][])>()
        {(new double[]{ 0.05301208792514899, 0.053825572057669926, 0.052065045951539365 }, "CenterNodeGradients",-1, new double[3][])};

        static double[] monitoredGPCoordsVelocity = new double[] { 0.09, 0.09, 0.09 };
        static double[] monitoredGPCoordsDivVelocity = new double[] { 0.09, 0.09, 0.09 };

        #endregion

        #region Darcy model

        //Darcy model properties

        // original case
        //static double Sv = 7e+3; // 1/(m)
        //static double Lp = 2.7e-12; // m/(KPa sec)
        //static double LplSvl_tumor = 0; // 1/(KPa sec)
        //static double LplSvl_host = 3.75e-1; // 1/(KPa sec)
        //static double pv = 4; // kPa
        //static double pl = 0d; // KPa
        //static double k_th = 7.52e-10; // m2/(KPa sec)
        //static double k_th_tumor = 7.52e-11; // m2/(KPa sec)
        //static double k_th_host = 7.52e-13; // m2/(KPa sec)

        // simplified case
        static double Sv = 0; // 1/(m)
        static double Lp = 0; // m/(KPa sec)
        static double LplSvl_tumor = 0; // 1/(KPa sec)
        static double LplSvl_host = 0; // 1/(KPa sec)
        static double pv = 0; // kPa
        static double pl = 0d; // KPa
        static double k_th_tumor = 7.52e-6; // m2/(KPa sec)
        static double k_th_host = 7.52e-6; // m2/(KPa sec)


        #endregion

        #region Darcy BCs

        private static ConvectionDiffusionDof[] constrainedDofType = new ConvectionDiffusionDof[1] { ConvectionDiffusionDof.UnknownVariable };
        private static double[] boundaryValue = new double[1] { 0d };
        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])> pressureDirichletBC =
            new List<(BC, ConvectionDiffusionDof[], double[][], double[])>()
            {
                (BC.BottomDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
                (BC.TopDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0.1}}, boundaryValue),
                (BC.LeftDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
                (BC.RightDirichlet, constrainedDofType, new double[1][]{new double[3] {0.1,0,0}}, boundaryValue),
                (BC.FrontDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
                (BC.BackDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0.1,0}}, boundaryValue),

            };

        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])>
            pressureNeumannBC = new List<(BC, ConvectionDiffusionDof[], double[][], double[])>();

        private static double boundaryValueAllBoundaries = 0;
        #endregion

        #region Darcy Initial condition values
        // Data 1:RegionType, 2:Bcstype, 3: Region CaracteristicCoords Id, 4: Bc value
        static List<(int, int, double[][], double[])> eq78InitialConditionsList = new List<(int, int, double[][], double[])>() { (0, 0, new double[3][], new double[3]) };

        private static double initialCondition = 0; // TODO Orestis delete when obsolete

        #endregion

        #region Darcy logs

        #region ToDo Orestis log task 3 
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: logged Doftype, 4: found node Id, 5: results 
        /// </summary>
        static List<(double[], string, StructuralDof, int, double[])> nodePressureLogs = new List<(double[], string, StructuralDof, int, double[])>()
        {(new double[]{ 0.04930793848882013,0.04994681648346263,0.075 }, "CornerNodeTranslationZ.txt",StructuralDof.TranslationZ,-1, new double[0])};

        static double[] pressureMonitorNodeCoords = new double[] { 0.055, 0.0559, 0.05 };
        private static int pressureMonitorID;
        static ConvectionDiffusionDof pressureMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        #endregion

        #region ToDo Orestis log task 4
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: found element Id, 4: results 
        /// </summary>
        static List<(double[], string, int, double[][])> gpPressureGraadientLogs = new List<(double[], string, int, double[][])>()
        {(new double[]{ 0.05301208792514899, 0.053825572057669926, 0.052065045951539365 }, "CenterNodePressureGradients",-1, new double[3][])};

        //static double[] monitoredGPcoordsPresGradient = new double[] { 0.004086132769345323, 0.006191771651571988, 0.09369050147682026 };
        static double[] monitoredGPcoordsPresGradient = new double[] { 0.055, 0.0559, 0.05 };

        #endregion
        #endregion

        #region Cancer Cell Density (TCell) model

        private const double dummySolidVelovity = -5;

        /// <summary>
        /// Growth rate parameter 1[mol/(m3)]
        /// </summary>
        private const double K1 = 1.74E-6; // [1/s]

        /// <summary>
        /// Growth rate parameter 2[mol/(m3)]
        /// </summary>
        private const double K2 = 8.3E-3; // [mol/(m3)]

        /// <summary>
        /// Oxygen concentration (Dependent Variable) [mol/m3]
        /// </summary>
        private const double Cox = 0.2; // [mol / m3]

        #endregion

        #region Cancer Cell Density (TCell) Boundary Conditions


        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])> tCellDirichletBC =
            new List<(BC, ConvectionDiffusionDof[], double[][], double[])>()
        {(BC.TopRightBackDiriclet, constrainedDofType, new double[2][]{new double[3] {0,0,0},new double[3] {0.1,0.1,0.1}}, new double[]{500d}),};

        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])>
                tCellNeumannBC = new List<(BC, ConvectionDiffusionDof[], double[][], double[])>();

        #endregion

        #region Cancer Cell Density (TCell)  Initial condition

        private double initialTCellDensity = 0d;

        #endregion

        #region Cancer Cell Density (TCell)  logs

        static List<(double[], string, StructuralDof, int, double[])> nodeTCellLogs = new List<(double[], string, StructuralDof, int, double[])>()
            {(new double[]{ 0.04930793848882013,0.04994681648346263,0.075 }, "CornerNodeTranslationZ.txt",StructuralDof.TranslationZ,-1, new double[0])};

        //static double[] tCellMonitorNodeCoords = new double[] { 0.055, 0.0559, 0.07366 };
        static double[] tCellMonitorNodeCoords = { 0.0, 0.0, 0.09 };

        private static int tCellMonitorID;

        static ConvectionDiffusionDof tCellMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        #endregion

        public UPTSolution()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        [Theory]

        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]

        public void MonophasicEquationModel(string fileName)
        {
            ContinuumElement3DGrowth.dT = timeStep;

            //Read geometry
            var comsolReader = new ComsolMeshReader(fileName);

            // initialize Shared quantities of Coupled model
            Dictionary<int, double> lambda = new Dictionary<int, double>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                lambda.Add(elem.Key, lambda0);
            }

            Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints =
                new Dictionary<int, double[][]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var gpTensorDiv = new double[nGaussPoints][];
                for (int i1 = 0; i1 < nGaussPoints; i1++)
                {
                    gpTensorDiv[i1] = new double[] { initial_dp_dx, initial_dp_dy, initial_dp_dz };
                }

                pressureTensorDivergenceAtElementGaussPoints.Add(elem.Key, gpTensorDiv);
            }

            Dictionary<int, double[]> velocityDivergenceAtElementGaussPoints =
                new Dictionary<int, double[]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var velocityDiv = new double[nGaussPoints];
                for (int i1 = 0; i1 < nGaussPoints; i1++)
                {
                    velocityDiv[i1] = velocityDivInitialVal;
                }

                velocityDivergenceAtElementGaussPoints.Add(elem.Key, velocityDiv);
            }

            Dictionary<int, double[][]> solidVelocity =
                new Dictionary<int, double[][]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var velocity = new double[nGaussPoints][];
                for (int i1 = 0; i1 < nGaussPoints; i1++)
                {
                    velocity[i1] = new double[] { 0d, 0d, 0d };
                }
                solidVelocity.Add(elem.Key, velocity);
            }

            Dictionary<int, double> dummyFieldCOx =
                new Dictionary<int, double>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                dummyFieldCOx.Add(elem.Key, Cox);
            }


            miNormal = miTumor; // TODO : remove this from here
            kappaNormal = kappaTumor;

            #region loggin (defined before model builder creation to give them nodes)

            structuralMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, structuralMonitorNodeCoords, 1e-2);
            pressureMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, pressureMonitorNodeCoords, 1e-2);
            tCellMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, tCellMonitorNodeCoords, 1e-2);

            var p_i = new double[(int)(totalTime / timeStep)];

            double[] structuralResultsX = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsY = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsZ = new double[(int)(totalTime / timeStep)];
            var displacements = new List<double[]>();
            displacements.Add(structuralResultsX);
            displacements.Add(structuralResultsY);
            displacements.Add(structuralResultsZ);
            
            double[] solidVelocityResultsX = new double[(int)(totalTime / timeStep)];
            double[] solidVelocityResultsY = new double[(int)(totalTime / timeStep)];
            double[] solidVelocityResultsZ = new double[(int)(totalTime / timeStep)];
            var solidVelocities = new List<double[]>();
            solidVelocities.Add(solidVelocityResultsX);
            solidVelocities.Add(solidVelocityResultsY);
            solidVelocities.Add(solidVelocityResultsZ);
            
            
            double[] gp_dut_dx_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dvt_dy_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dwt_dz_OverTime = new double[(int)(totalTime / timeStep)];
            var divVelocity = new List<double[]>();
            divVelocity.Add(gp_dut_dx_OverTime);
            divVelocity.Add(gp_dvt_dy_OverTime);
            divVelocity.Add(gp_dwt_dz_OverTime);


            double[] gp_div_v_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dx_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dy_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dz_Overtime = new double[(int)(totalTime / timeStep)];
            var dp_dxi = new List<double[]>();
            dp_dxi.Add(gp_dP_dx_OverTime);
            dp_dxi.Add(gp_dP_dy_OverTime);
            dp_dxi.Add(gp_dP_dz_Overtime);


            double[] tCell = new double[(int)(totalTime / timeStep)];

            int monitoredGPDivVelocity_elemID = -1; 
            int monitoredGPpressureGrad_elemID = -1; 
            int monitoredGPVelocity_elemID = -1;

            //TODO Orestis: implement here one for loop for each Gp Type requested Output
            #endregion

            //Create Model For Pressure
            var pressureModel = new Eq78ModelProviderForStaggeredSolutionex7ref(comsolReader, k_th_tumor, k_th_host, Lp, Sv, pv,
                LplSvl_tumor, LplSvl_host, pl, velocityDivergenceAtElementGaussPoints, pressureMonitorID, pressureMonitorDOF, pressureDirichletBC, pressureNeumannBC);

            //Create Model For Structural
            var structuralModel = new Eq9ModelProviderForStaggeredSolutionEx7Ref(comsolReader, Sc, miNormal, kappaNormal, miTumor,
                kappaTumor, density, timeStep, totalTime, lambda, pressureTensorDivergenceAtElementGaussPoints,
                structuralMonitorID, structuralMonitorDOF, structuralNeumannBC, structuralDirichletBC);


            //Create Model For TCell
            var tCellModel = new TCellModelProvider(K1, K2, dummyFieldCOx, solidVelocity, comsolReader, tCellMonitorDOF, tCellMonitorID, tCellDirichletBC, tCellNeumannBC, initialTCellDensity);
            

            var equationModel = new CoupledBiphasicTCellModelProvider(pressureModel, structuralModel, tCellModel, comsolReader, lambda,
                pressureTensorDivergenceAtElementGaussPoints, velocityDivergenceAtElementGaussPoints, solidVelocity,timeStep,
                totalTime, incrementsPertimeStep);

            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers,
                equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 200, tolerance: 1E-5);
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                equationModel.CurrentTimeStep = currentTimeStep;
                equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
                staggeredAnalyzer.SolveCurrentStep();

                #region logging
                
                monitoredGPDivVelocity_elemID = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], monitoredGPCoordsDivVelocity, 1e-1); //Todo Orestis delete these commands1
                monitoredGPpressureGrad_elemID = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], monitoredGPcoordsPresGradient, 1e-1);
                monitoredGPVelocity_elemID = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], monitoredGPCoordsVelocity, 1e-1); //Todo Orestis delete these commands1

                //nodal logs
                p_i[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[0].GetNode(pressureMonitorID), pressureMonitorDOF];
                //p_i[currentTimeStep] = 0d;
                //structuralResultsX[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationX];
                //structuralResultsY[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationY];
                //structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationZ];
                structuralResultsX[currentTimeStep] = 0d;
                structuralResultsY[currentTimeStep] = 0d;
                structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), structuralMonitorDOF];

                //gp (element) logs
                solidVelocityResultsX[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocity[0][0];
                solidVelocityResultsY[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocity[0][1];
                solidVelocityResultsZ[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocity[0][2];
                
                gp_dP_dx_OverTime[currentTimeStep] = ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGPpressureGrad_elemID]).xcoeff_OverTimeAtGp1[0];
                gp_dP_dy_OverTime[currentTimeStep] = ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGPpressureGrad_elemID]).ycoeff_OverTimeAtGp1[0];
                gp_dP_dz_Overtime[currentTimeStep] = ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGPpressureGrad_elemID]).zcoeff_OverTimeAtGp1[0];
                
                gp_dut_dx_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPDivVelocity_elemID]).velocityDivergence_term1[0];
                gp_dvt_dy_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPDivVelocity_elemID]).velocityDivergence_term2[0];
                gp_dwt_dz_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPDivVelocity_elemID]).velocityDivergence_term3[0];
                
                gp_div_v_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPDivVelocity_elemID]).velocityDivergence[0];
                

                tCell[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[2].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[2].GetNode(tCellMonitorID), tCellMonitorDOF];

                //model maximus (DO NOT ERASE)
                //modelMaxVelDivOverTime[currentTimeStep] = velocityDivergenceAtElementGaussPoints.Select(x => Math.Abs(x.Value[0])).ToArray().Max();
                //modelMax_dP_dxOverTime[currentTimeStep] = pressureTensorDivergenceAtElementGaussPoints.Select(x => Math.Abs(x.Value[0][0])).ToArray().Max();


                /*if (Solution.ContainsKey(currentTimeStep))
                {
                    Solution[currentTimeStep] = allValues;
                    Console.WriteLine($"Time step: {timeStep}");
                    Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[timeStep])}");
                }
                else
                {
                    Solution.Add(currentTimeStep, allValues);
                }*/

                #endregion
                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    (equationModel.ParentAnalyzers[j] as NewmarkDynamicAnalyzer).AdvanceStep();

                }
                
                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
                    equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
                }

                equationModel.SaveStateFromElements();

                //Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[currentTimeStep])}");
            }

            Assert.True(ResultChecker.CheckResults(structuralResultsZ, expectedDisplacments(), 1e-3));
            
            Assert.True(ResultChecker.CheckResults(solidVelocityResultsX, expectedSolidVelocityX(), 1e-3));
            Assert.True(ResultChecker.CheckResults(solidVelocityResultsY, expectedSolidVelocityY(), 1e-3));
            Assert.True(ResultChecker.CheckResults(solidVelocityResultsZ, expectedSolidVelocityZ(), 1e-3));
            
            Assert.True(ResultChecker.CheckResults(p_i, expectedPressurevalues(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dut_dx_OverTime, expected_dutdx_values(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dvt_dy_OverTime, expected_dvtdy_values(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dwt_dz_OverTime, expected_dwtdz_values(), 1e-3));
            //Assert.True(ResultChecker.CheckResults(gp_div_v_OverTime, expected_div_vs_values(), 1e-1));
            Assert.True(ResultChecker.CheckResults(gp_dP_dx_OverTime, expected_dpdx_values(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dP_dy_OverTime, expected_dpdy_values(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dP_dz_Overtime, expected_dpdz_values(), 1e-3));
            //Assert.True(ResultChecker.CheckResults(tCell, expected_Tc_values(), 1e-3));





            double pr = 100;
            double Fval = 100;
            double Fval_e = 0;
            //eq9model.load_value;

            //pressure_{ pr}_F_{Fval}_e{Fval_e}_LOGGEDval_ PAth name do not erase this exampleNo_{exNo}_caseNo_{caseNo}_
            /*var writer = new MGroup.LinearAlgebra.Output.Array1DWriter();
            var patSelection = 0;
            var outputPath = patSelection == 0 ? "../../../StaggeredSolutionPresDynamex7ref/results1/" : $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\";
            writer.WriteToFile(structuralResultsX, outputPath + $@"u_{1}_.txt");
            writer.WriteToFile(structuralResultsY, outputPath + $@"u_{2}_.txt");
            writer.WriteToFile(structuralResultsZ, outputPath + $@"u_{3}_.txt");
            writer.WriteToFile(gp_dut_dx_OverTime, outputPath + $@"gp_dut_dx_OverTime.txt");
            writer.WriteToFile(gp_dvt_dy_OverTime, outputPath + $@"gp_dvt_dy_OverTime.txt");
            writer.WriteToFile(gp_dwt_dz_OverTime, outputPath + $@"gp_dwt_dz_OverTime.txt");
            writer.WriteToFile(gp_div_v_OverTime, outputPath +  $@"gp_div_v_OverTime_.txt");



            writer.WriteToFile(p_i, outputPath + $@"p_i.txt");
            writer.WriteToFile(gp_dP_dx_OverTime, outputPath + $@"gp_dP_dx_OverTime.txt");
            writer.WriteToFile(gp_dP_dy_OverTime, outputPath + $@"gp_dP_dy_OverTime.txt");
            writer.WriteToFile(gp_dP_dz_Overtime, outputPath + $@"gp_dP_dz_Overtime.txt");*/


            //var path = outputPath+"dp_dxi_mslv.csv";
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(dp_dxi), "../../../StaggeredBiphasicTCell/dp_dxi_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(displacements), "../../../StaggeredBiphasicTCell/displacements_nodes_mslv.csv");
            CSVExporter.ExportVectorToCSV(p_i, "../../../StaggeredBiphasicTCell/pi_nodes_mslv.csv");
            CSVExporter.ExportVectorToCSV(tCell, "../../../StaggeredBiphasicTCell/tCell_nodes_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(divVelocity), "../../../StaggeredBiphasicTCell/dut_dxi_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(solidVelocities), "../../../StaggeredBiphasicTCell/ut_GP_mslv.csv");

        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //Node : 0.0, 0.0, 0.1
        public static double[] expectedDisplacments()
        {
            return new double[] {
                -3.66177E-09,
                -1.82506E-08,
                -4.72353E-08,
                -9.03487E-08,
                -1.47347E-07,
                -2.18019E-07,
                -3.02169E-07,
                -3.99613E-07,
                -5.10167E-07,
                -6.33648E-07,
            }; 
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //GP : 0.09, 0.09, 0.09
        private static double[] expectedSolidVelocityX()
        {
            return new double[] {
                -2.57948E-08,
                -1.5538E-07 ,
                -4.79338E-07,
                -1.00336E-06,
                -1.65214E-06,
                -2.38495E-06,
                -3.20718E-06,
                -4.12942E-06,
                -5.16614E-06,
                -6.33261E-06
            };
        }
        
        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //GP : 0.09, 0.09, 0.09
        private static double[] expectedSolidVelocityY()
        {
            return new double[] {
                -5.30026E-08,
                -3.23281E-07,
                -1.00454E-06,
                -2.11387E-06,
                -3.49449E-06,
                -5.04341E-06,
                -6.74082E-06,
                -8.57872E-06,
                -1.05588E-05,
                -1.26869E-05
            };
        }
        
        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //GP : 0.09, 0.09, 0.09
        private static double[] expectedSolidVelocityZ()
        {
            return new double[] {
                -3.31285E-05,
                -0.000132685,
                -0.000265916,
                -0.000399804,
                -0.000534109,
                -0.000668626,
                -0.000803237,
                -0.000937851,
                -0.001072391,
                -0.001206793,
            };
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //GP : 0.05, 0.05, 0.05
        public static double[] expectedPressurevalues()
        {
            return new double[] {
                -2.62803E-05,
                -7.73618E-05,
                -0.000142106,
                -0.000178119,
                -0.000222339,
                -0.000265597,
                -0.000295448,
                -0.000282738,
                -0.00025486 ,
                -0.000182668
            };
        }
        
        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.09, 0.09, 0.09
        public static double[] expected_dutdx_values()
        {
            return new double[] {
                3.09537E-06,
                1.86456E-05,
                5.75206E-05,
                0.000120403,
                0.000198256,
                0.000286194,
                0.000384862,
                0.000495531,
                0.000619937,
                0.000759913,
            };
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.09, 0.09, 0.09
        public static double[] expected_dvtdy_values()
        {
            return new double[] {
                6.36031E-06,
                3.87937E-05,
                0.000120545,
                0.000253664,
                0.000419339,
                0.000605209,
                0.000808898,
                0.001029447,
                0.001267052,
                0.001522424,
            };
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.09, 0.09, 0.09
        public static double[] expected_dwtdz_values()
        {
            return new double[] {
                -0.006214244,
                -0.024834619,
                -0.049594341,
                -0.07424962 ,
                -0.098808016,
                -0.123276063,
                -0.147651799,
                -0.17192928 ,
                -0.196099345,
                -0.220151411,
            };
        }

        public static double[] expected_div_vs_values()
        {
            return new double[10];
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.05, 0.05, 0.05
        public static double[] expected_dpdx_values()
        {
            return new double[] {
                -0.000510087,
                -0.000582082,
                0.000128969 ,
                0.005425463 ,
                0.011910756 ,
                0.020264921 ,
                0.028152223 ,
                0.033841904 ,
                0.038328132 ,
                0.039359541 ,
            };
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.05, 0.05, 0.05
        public static double[] expected_dpdy_values()
        {
            return new double[]
            {
                -0.001559998,
                -0.003938722,
                -0.006777791,
                -0.007825571,
                -0.010430786,
                -0.017679464,
                -0.026326288,
                -0.039003849,
                -0.052556255,
                -0.068290536,
            };
        }
        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.05, 0.05, 0.05
        public static double[] expected_dpdz_values()
        {
            return new double[]
            {
                -0.00085234,
                -0.001205138,
                -0.001344107,
                1.11914E-05,
                -0.000479235,
                -0.007219547,
                -0.016283798,
                -0.031968813,
                -0.049267171,
                -0.07062152,
            };
        }

        public static double[] expected_Tc_values()
        {
            return new double[] {
                8.97E-03,
                3.19E-02,
                7.16E-02,
                1.30E-01,
                2.09E-01,
                3.08E-01,
                4.29E-01,
                5.72E-01,
                7.39E-01,
                9.29E-01,
            };
        }







    }
}
