/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Sun Apr 5 17:58:46 2015
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QStatusBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "glwidget.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionExit;
    QAction *actionReset_Camera;
    QAction *actionTake_Screenshot;
    QAction *actionExport_OBJ;
    QAction *actionImport_OBJ;
    QAction *actionImport_Metric;
    QAction *actionAdd_Noise;
    QAction *actionSet_No_Target_Metric;
    QAction *actionSet_Negative_K_Target_Metric;
    QAction *actionMinimize_with_Newton;
    QAction *actionSymmetrize;
    QAction *actionEigenvalues;
    QAction *actionMake_Cone;
    QAction *actionMake_Flat_Cone;
    QAction *actionSet_Current_Lengths_as_Intrinsic;
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout;
    GLWidget *GLwidget;
    QWidget *widget;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout_9;
    QGroupBox *renderingBox;
    QCheckBox *smoothShadeCheckBox;
    QCheckBox *wireframeCheckBox;
    QLabel *outputLabel;
    QLineEdit *outputEdit;
    QGroupBox *algorithmsBox;
    QPushButton *crushButton;
    QGroupBox *parameterBox;
    QGroupBox *physicalBox;
    QLabel *youngsModulusLabel;
    QLineEdit *youngsModulusEdit;
    QLabel *poissonRatioLabel;
    QLineEdit *poissonRatioEdit;
    QLabel *thicknessLabel;
    QLineEdit *thicknessEdit;
    QLabel *densityLabel;
    QLineEdit *densityEdit;
    QLabel *scaleLabel;
    QLineEdit *scaleEdit;
    QGroupBox *solverBox;
    QLabel *eulerItersLabel;
    QLabel *eulerTimestepLabel;
    QLineEdit *eulerItersEdit;
    QLineEdit *eulerTimestepEdit;
    QLabel *dampingCoeffLabel;
    QLineEdit *dampingCoeffEdit;
    QGroupBox *problemBox;
    QLabel *growthAmountLabel;
    QLabel *baseProbabilityLabel;
    QLabel *maxStrainLabel;
    QLineEdit *growthAmountEdit;
    QLineEdit *maxStrainEdit;
    QLineEdit *baseProbabilityEdit;
    QMenuBar *menubar;
    QMenu *menuVViewer;
    QMenu *menuView;
    QMenu *menuActions;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1600, 1200);
        MainWindow->setBaseSize(QSize(0, 0));
        actionExit = new QAction(MainWindow);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionReset_Camera = new QAction(MainWindow);
        actionReset_Camera->setObjectName(QString::fromUtf8("actionReset_Camera"));
        actionTake_Screenshot = new QAction(MainWindow);
        actionTake_Screenshot->setObjectName(QString::fromUtf8("actionTake_Screenshot"));
        actionExport_OBJ = new QAction(MainWindow);
        actionExport_OBJ->setObjectName(QString::fromUtf8("actionExport_OBJ"));
        actionImport_OBJ = new QAction(MainWindow);
        actionImport_OBJ->setObjectName(QString::fromUtf8("actionImport_OBJ"));
        actionImport_Metric = new QAction(MainWindow);
        actionImport_Metric->setObjectName(QString::fromUtf8("actionImport_Metric"));
        actionAdd_Noise = new QAction(MainWindow);
        actionAdd_Noise->setObjectName(QString::fromUtf8("actionAdd_Noise"));
        actionSet_No_Target_Metric = new QAction(MainWindow);
        actionSet_No_Target_Metric->setObjectName(QString::fromUtf8("actionSet_No_Target_Metric"));
        actionSet_Negative_K_Target_Metric = new QAction(MainWindow);
        actionSet_Negative_K_Target_Metric->setObjectName(QString::fromUtf8("actionSet_Negative_K_Target_Metric"));
        actionMinimize_with_Newton = new QAction(MainWindow);
        actionMinimize_with_Newton->setObjectName(QString::fromUtf8("actionMinimize_with_Newton"));
        actionSymmetrize = new QAction(MainWindow);
        actionSymmetrize->setObjectName(QString::fromUtf8("actionSymmetrize"));
        actionEigenvalues = new QAction(MainWindow);
        actionEigenvalues->setObjectName(QString::fromUtf8("actionEigenvalues"));
        actionMake_Cone = new QAction(MainWindow);
        actionMake_Cone->setObjectName(QString::fromUtf8("actionMake_Cone"));
        actionMake_Flat_Cone = new QAction(MainWindow);
        actionMake_Flat_Cone->setObjectName(QString::fromUtf8("actionMake_Flat_Cone"));
        actionSet_Current_Lengths_as_Intrinsic = new QAction(MainWindow);
        actionSet_Current_Lengths_as_Intrinsic->setObjectName(QString::fromUtf8("actionSet_Current_Lengths_as_Intrinsic"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        horizontalLayout = new QHBoxLayout(centralwidget);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        GLwidget = new GLWidget(centralwidget);
        GLwidget->setObjectName(QString::fromUtf8("GLwidget"));
        GLwidget->setMinimumSize(QSize(100, 100));

        horizontalLayout->addWidget(GLwidget);

        widget = new QWidget(centralwidget);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setMinimumSize(QSize(310, 310));
        widget->setMaximumSize(QSize(190, 16777215));
        verticalLayoutWidget = new QWidget(widget);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(50, 10, 220, 704));
        verticalLayout_9 = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout_9->setObjectName(QString::fromUtf8("verticalLayout_9"));
        verticalLayout_9->setContentsMargins(0, 0, 0, 0);
        renderingBox = new QGroupBox(verticalLayoutWidget);
        renderingBox->setObjectName(QString::fromUtf8("renderingBox"));
        renderingBox->setMinimumSize(QSize(0, 40));
        renderingBox->setMaximumSize(QSize(16777215, 16777215));
        smoothShadeCheckBox = new QCheckBox(renderingBox);
        smoothShadeCheckBox->setObjectName(QString::fromUtf8("smoothShadeCheckBox"));
        smoothShadeCheckBox->setGeometry(QRect(10, 20, 131, 22));
        smoothShadeCheckBox->setChecked(true);
        wireframeCheckBox = new QCheckBox(renderingBox);
        wireframeCheckBox->setObjectName(QString::fromUtf8("wireframeCheckBox"));
        wireframeCheckBox->setGeometry(QRect(10, 40, 151, 22));
        wireframeCheckBox->setChecked(true);
        outputLabel = new QLabel(renderingBox);
        outputLabel->setObjectName(QString::fromUtf8("outputLabel"));
        outputLabel->setGeometry(QRect(10, 70, 121, 17));
        outputEdit = new QLineEdit(renderingBox);
        outputEdit->setObjectName(QString::fromUtf8("outputEdit"));
        outputEdit->setGeometry(QRect(140, 70, 61, 20));

        verticalLayout_9->addWidget(renderingBox);

        algorithmsBox = new QGroupBox(verticalLayoutWidget);
        algorithmsBox->setObjectName(QString::fromUtf8("algorithmsBox"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(algorithmsBox->sizePolicy().hasHeightForWidth());
        algorithmsBox->setSizePolicy(sizePolicy);
        algorithmsBox->setMinimumSize(QSize(0, 60));
        algorithmsBox->setSizeIncrement(QSize(0, 0));
        algorithmsBox->setBaseSize(QSize(0, 0));
        crushButton = new QPushButton(algorithmsBox);
        crushButton->setObjectName(QString::fromUtf8("crushButton"));
        crushButton->setGeometry(QRect(60, 60, 99, 27));

        verticalLayout_9->addWidget(algorithmsBox);

        parameterBox = new QGroupBox(verticalLayoutWidget);
        parameterBox->setObjectName(QString::fromUtf8("parameterBox"));
        parameterBox->setMinimumSize(QSize(0, 500));
        physicalBox = new QGroupBox(parameterBox);
        physicalBox->setObjectName(QString::fromUtf8("physicalBox"));
        physicalBox->setGeometry(QRect(10, 30, 191, 171));
        physicalBox->setAutoFillBackground(false);
        physicalBox->setFlat(false);
        youngsModulusLabel = new QLabel(physicalBox);
        youngsModulusLabel->setObjectName(QString::fromUtf8("youngsModulusLabel"));
        youngsModulusLabel->setGeometry(QRect(0, 60, 121, 17));
        youngsModulusEdit = new QLineEdit(physicalBox);
        youngsModulusEdit->setObjectName(QString::fromUtf8("youngsModulusEdit"));
        youngsModulusEdit->setGeometry(QRect(130, 60, 61, 20));
        poissonRatioLabel = new QLabel(physicalBox);
        poissonRatioLabel->setObjectName(QString::fromUtf8("poissonRatioLabel"));
        poissonRatioLabel->setGeometry(QRect(0, 90, 121, 17));
        poissonRatioEdit = new QLineEdit(physicalBox);
        poissonRatioEdit->setObjectName(QString::fromUtf8("poissonRatioEdit"));
        poissonRatioEdit->setGeometry(QRect(130, 90, 61, 20));
        thicknessLabel = new QLabel(physicalBox);
        thicknessLabel->setObjectName(QString::fromUtf8("thicknessLabel"));
        thicknessLabel->setGeometry(QRect(0, 120, 121, 17));
        thicknessEdit = new QLineEdit(physicalBox);
        thicknessEdit->setObjectName(QString::fromUtf8("thicknessEdit"));
        thicknessEdit->setGeometry(QRect(130, 120, 61, 20));
        densityLabel = new QLabel(physicalBox);
        densityLabel->setObjectName(QString::fromUtf8("densityLabel"));
        densityLabel->setGeometry(QRect(0, 150, 67, 17));
        densityEdit = new QLineEdit(physicalBox);
        densityEdit->setObjectName(QString::fromUtf8("densityEdit"));
        densityEdit->setGeometry(QRect(130, 150, 61, 20));
        scaleLabel = new QLabel(physicalBox);
        scaleLabel->setObjectName(QString::fromUtf8("scaleLabel"));
        scaleLabel->setGeometry(QRect(0, 30, 121, 17));
        scaleEdit = new QLineEdit(physicalBox);
        scaleEdit->setObjectName(QString::fromUtf8("scaleEdit"));
        scaleEdit->setGeometry(QRect(130, 30, 61, 20));
        solverBox = new QGroupBox(parameterBox);
        solverBox->setObjectName(QString::fromUtf8("solverBox"));
        solverBox->setGeometry(QRect(10, 210, 191, 121));
        eulerItersLabel = new QLabel(solverBox);
        eulerItersLabel->setObjectName(QString::fromUtf8("eulerItersLabel"));
        eulerItersLabel->setGeometry(QRect(0, 30, 111, 17));
        eulerTimestepLabel = new QLabel(solverBox);
        eulerTimestepLabel->setObjectName(QString::fromUtf8("eulerTimestepLabel"));
        eulerTimestepLabel->setGeometry(QRect(0, 60, 101, 17));
        eulerItersEdit = new QLineEdit(solverBox);
        eulerItersEdit->setObjectName(QString::fromUtf8("eulerItersEdit"));
        eulerItersEdit->setGeometry(QRect(130, 30, 61, 20));
        eulerTimestepEdit = new QLineEdit(solverBox);
        eulerTimestepEdit->setObjectName(QString::fromUtf8("eulerTimestepEdit"));
        eulerTimestepEdit->setGeometry(QRect(130, 60, 61, 20));
        dampingCoeffLabel = new QLabel(solverBox);
        dampingCoeffLabel->setObjectName(QString::fromUtf8("dampingCoeffLabel"));
        dampingCoeffLabel->setGeometry(QRect(0, 90, 111, 17));
        dampingCoeffEdit = new QLineEdit(solverBox);
        dampingCoeffEdit->setObjectName(QString::fromUtf8("dampingCoeffEdit"));
        dampingCoeffEdit->setGeometry(QRect(130, 90, 61, 20));
        problemBox = new QGroupBox(parameterBox);
        problemBox->setObjectName(QString::fromUtf8("problemBox"));
        problemBox->setGeometry(QRect(10, 330, 201, 171));
        growthAmountLabel = new QLabel(problemBox);
        growthAmountLabel->setObjectName(QString::fromUtf8("growthAmountLabel"));
        growthAmountLabel->setGeometry(QRect(0, 60, 111, 17));
        baseProbabilityLabel = new QLabel(problemBox);
        baseProbabilityLabel->setObjectName(QString::fromUtf8("baseProbabilityLabel"));
        baseProbabilityLabel->setGeometry(QRect(0, 120, 121, 17));
        maxStrainLabel = new QLabel(problemBox);
        maxStrainLabel->setObjectName(QString::fromUtf8("maxStrainLabel"));
        maxStrainLabel->setGeometry(QRect(0, 90, 111, 17));
        growthAmountEdit = new QLineEdit(problemBox);
        growthAmountEdit->setObjectName(QString::fromUtf8("growthAmountEdit"));
        growthAmountEdit->setGeometry(QRect(130, 60, 61, 20));
        maxStrainEdit = new QLineEdit(problemBox);
        maxStrainEdit->setObjectName(QString::fromUtf8("maxStrainEdit"));
        maxStrainEdit->setGeometry(QRect(130, 90, 61, 20));
        baseProbabilityEdit = new QLineEdit(problemBox);
        baseProbabilityEdit->setObjectName(QString::fromUtf8("baseProbabilityEdit"));
        baseProbabilityEdit->setGeometry(QRect(130, 120, 61, 20));

        verticalLayout_9->addWidget(parameterBox);


        horizontalLayout->addWidget(widget);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 1600, 25));
        menuVViewer = new QMenu(menubar);
        menuVViewer->setObjectName(QString::fromUtf8("menuVViewer"));
        menuView = new QMenu(menubar);
        menuView->setObjectName(QString::fromUtf8("menuView"));
        menuActions = new QMenu(menubar);
        menuActions->setObjectName(QString::fromUtf8("menuActions"));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);

        menubar->addAction(menuVViewer->menuAction());
        menubar->addAction(menuView->menuAction());
        menubar->addAction(menuActions->menuAction());
        menuVViewer->addAction(actionImport_OBJ);
        menuVViewer->addAction(actionImport_Metric);
        menuVViewer->addAction(actionExit);
        menuView->addAction(actionReset_Camera);
        menuView->addAction(actionTake_Screenshot);
        menuView->addAction(actionExport_OBJ);
        menuActions->addSeparator();
        menuActions->addAction(actionMake_Cone);
        menuActions->addAction(actionMake_Flat_Cone);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("MainWindow", "Exit", 0, QApplication::UnicodeUTF8));
        actionReset_Camera->setText(QApplication::translate("MainWindow", "Reset Camera", 0, QApplication::UnicodeUTF8));
        actionTake_Screenshot->setText(QApplication::translate("MainWindow", "Take Screenshot", 0, QApplication::UnicodeUTF8));
        actionExport_OBJ->setText(QApplication::translate("MainWindow", "Export 3D Mesh", 0, QApplication::UnicodeUTF8));
        actionImport_OBJ->setText(QApplication::translate("MainWindow", "Import OBJ", 0, QApplication::UnicodeUTF8));
        actionImport_Metric->setText(QApplication::translate("MainWindow", "Import Metric", 0, QApplication::UnicodeUTF8));
        actionAdd_Noise->setText(QApplication::translate("MainWindow", "Add Noise", 0, QApplication::UnicodeUTF8));
        actionSet_No_Target_Metric->setText(QApplication::translate("MainWindow", "Set No Target Metric", 0, QApplication::UnicodeUTF8));
        actionSet_Negative_K_Target_Metric->setText(QApplication::translate("MainWindow", "Set Negative K Target Metric", 0, QApplication::UnicodeUTF8));
        actionMinimize_with_Newton->setText(QApplication::translate("MainWindow", "Minimize with Newton", 0, QApplication::UnicodeUTF8));
        actionSymmetrize->setText(QApplication::translate("MainWindow", "Symmetrize", 0, QApplication::UnicodeUTF8));
        actionEigenvalues->setText(QApplication::translate("MainWindow", "Eigenvalues", 0, QApplication::UnicodeUTF8));
        actionMake_Cone->setText(QApplication::translate("MainWindow", "Make Cone", 0, QApplication::UnicodeUTF8));
        actionMake_Flat_Cone->setText(QApplication::translate("MainWindow", "Make Flat Cone", 0, QApplication::UnicodeUTF8));
        actionSet_Current_Lengths_as_Intrinsic->setText(QApplication::translate("MainWindow", "Set Current Lengths as Intrinsic", 0, QApplication::UnicodeUTF8));
        renderingBox->setTitle(QApplication::translate("MainWindow", "Rendering", 0, QApplication::UnicodeUTF8));
        smoothShadeCheckBox->setText(QApplication::translate("MainWindow", "Smooth Shade", 0, QApplication::UnicodeUTF8));
        wireframeCheckBox->setText(QApplication::translate("MainWindow", "Show Wireframe", 0, QApplication::UnicodeUTF8));
        outputLabel->setText(QApplication::translate("MainWindow", "Output folder:", 0, QApplication::UnicodeUTF8));
        algorithmsBox->setTitle(QApplication::translate("MainWindow", "Algorithms", 0, QApplication::UnicodeUTF8));
        crushButton->setText(QApplication::translate("MainWindow", "Crush", 0, QApplication::UnicodeUTF8));
        parameterBox->setTitle(QApplication::translate("MainWindow", "Parameters", 0, QApplication::UnicodeUTF8));
        physicalBox->setTitle(QApplication::translate("MainWindow", "Physical", 0, QApplication::UnicodeUTF8));
        youngsModulusLabel->setText(QApplication::translate("MainWindow", "Young's Modulus:", 0, QApplication::UnicodeUTF8));
        poissonRatioLabel->setText(QApplication::translate("MainWindow", "Poisson Ratio:", 0, QApplication::UnicodeUTF8));
        thicknessLabel->setText(QApplication::translate("MainWindow", "Thickness:", 0, QApplication::UnicodeUTF8));
        densityLabel->setText(QApplication::translate("MainWindow", "Density:", 0, QApplication::UnicodeUTF8));
        scaleLabel->setText(QApplication::translate("MainWindow", "Scale:", 0, QApplication::UnicodeUTF8));
        solverBox->setTitle(QApplication::translate("MainWindow", "Solver", 0, QApplication::UnicodeUTF8));
        eulerItersLabel->setText(QApplication::translate("MainWindow", "Euler Iterations", 0, QApplication::UnicodeUTF8));
        eulerTimestepLabel->setText(QApplication::translate("MainWindow", "Euler Timestep", 0, QApplication::UnicodeUTF8));
        dampingCoeffLabel->setText(QApplication::translate("MainWindow", "Damping Coeff", 0, QApplication::UnicodeUTF8));
        problemBox->setTitle(QApplication::translate("MainWindow", "Problem", 0, QApplication::UnicodeUTF8));
        growthAmountLabel->setText(QApplication::translate("MainWindow", "Growth Amount", 0, QApplication::UnicodeUTF8));
        baseProbabilityLabel->setText(QApplication::translate("MainWindow", "Base Probability", 0, QApplication::UnicodeUTF8));
        maxStrainLabel->setText(QApplication::translate("MainWindow", "Max Strain", 0, QApplication::UnicodeUTF8));
        menuVViewer->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuView->setTitle(QApplication::translate("MainWindow", "View", 0, QApplication::UnicodeUTF8));
        menuActions->setTitle(QApplication::translate("MainWindow", "Actions", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
