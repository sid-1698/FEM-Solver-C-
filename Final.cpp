/* ******************************************FEM SOLVER***************************************************** */
/* *************************************  S.SIDARTH MDM16B032 ********************************************** */

//Header Files
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <conio.h>
using namespace std;

#define N 10
#define K 30

float element_stiffness[6][6];
float global_stiffness[K][K];
float modified_stiffnesss[K][K];

float global_force[K][1];
float global_force_final[K][1];
float modified_force[K][1];

float global_displacement[K][2];
float solution[K][1];


//Main Menu
char menu(){
    char ch;
    system("cls");
    cout<<"*****************MENU****************"<<endl;
    cout<<"1. Add Nodes"<<endl;
    cout<<"2. Define Elements"<<endl;
    cout<<"3. View Element Connectivity Table"<<endl;
    cout<<"4. Add Forces"<<endl;
    cout<<"5. Define Boundary Conditions"<<endl;
    cout<<"6. View Solution"<<endl;
    cout<<"7. View Secondary Variable (Stress and Reaction Forces)"<<endl;
    cout<<"8. Quit"<<endl;
    cout<<"Enter your Choice (1-6) ";
    ch = getch();

    return ch;
}//End of menu

//Set Matrices to Zero
void settozero(float a[6][6],float b[K][K], float c[K][1],int key){
    int i,j;
    if(key == 1){
        for(i=0;i<6;i++)
            for(j=0;j<6;j++)
                a[i][j]=0;
    }
    else
    if(key == 2){
        for(i=0;i<K;i++)
            for(j=0;j<K;j++)
                b[i][j]=0;
    }
    else
    if(key == 3){
        for(i=0;i<K;i++)
            c[i][0]=0;
    }
}//End of settozero

//View Element Connectivity Table
void View(float elements[N][8],int e,int n){
    system("cls");
    cout<<"********************************ELEMENT CONNECTIVITY TABLE********************************"<<endl;
    cout<<setw(10)<<"ELEMENT|"<<setw(10)<<"NODE i|"<<setw(10)<<"NODE j|"<<setw(10)<<"LENGTH|"<<setw(15)<<"E(Pa)|";
    cout<<setw(15)<<"A(m^2)|"<<setw(15)<<"I(m^4)|"<<setw(10)<<"l(cos)|"<<setw(10)<<"m(sin)|"<<setw(10)<<"Alpha(/C)"<<endl;
    int i,j;
    for(i=0;i<115;i++){
        cout<<"-";
    }
    cout<<endl;
    for(i=0;i<e;i++){
        cout<<setw(9)<<i+1<<"|";
        for(j=0;j<8;j++){
            if(j<=2||j>=6)
                cout<<setw(9)<<setprecision(5)<<elements[i][j]<<"|";
            else
            if(j>2&&j<6)
                cout<<setw(14)<<setprecision(5)<<elements[i][j]<<"|";
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"Number of Elements : "<<e<<endl;
    cout<<"Number of Nodes : "<<n<<endl;
}//End of view

//View the Nodes with its coordinates
void View_nodes(int nodes[N][2],int n){
    system("cls");
    cout<<setw(10)<<"NODE NO|"<<setw(10)<<"X|"<<setw(9)<<"Y"<<endl;
    int i;
    for(i=0;i<30;i++){
        cout<<"-";
    }
    cout<<endl;
    for(i=0;i<n;i++){
        cout<<setw(9)<<setprecision(9)<<i+1<<"|"<<setw(9)<<setprecision(9)<<nodes[i][0]<<"|"<<setw(9)<<nodes[i][1]<<endl;
    }
}//End of View_nodes

//Check if the node already exists
int checknode(int nodes[N][2],int p){

    int x=nodes[p][0];
    int y=nodes[p][1];
    int i;
    for(i=0;i<p;i++){
        if(nodes[i][0]==x&&nodes[i][1]==y)
        return 0;
    }
    return 1;
}//End of checknode

//Add Nodes
void Nodes(int nodes[N][2],int *n){
    int i = *n;
    if(i==0){
        system("cls");
        cout<<"No nodes exist!"<<endl;
    }
    else{
        View_nodes(nodes,i);
    }
    char ch;
    cout<<"\n\n";
    cout<<"Press B to go back. Enter to Add a Node";
    ch = getch();
    if(ch == 'B'||ch == 'b'){
        return;
    }
    cout<<endl;
    cout<<"Enter the coordinates of Node "<<i+1<<endl;
    cout<<"Enter X coordinate :"; cin>>nodes[i][0];
    cout<<"Enter Y coordinate :"; cin>>nodes[i][1];
    if(checknode(nodes,i)==0){
        cout<<"Sorry! Node Point already exists.";
    }
    else{
        cout<<"Node Successfully Added.";
        *n=i+1;
    }
    cout<<"Press any key to go back.";
    getch();
}//End of Nodes

//Checks if element already exists
int checkelement(float elements[N][8],int p,int n){
    int i=elements[p][0];
    int j=elements[p][1];
    int k;
    if(i>n||j>n){
        return 0;
    }
    return 1;
}//End of checkelement

//Compute element length
float Length(int nodes[N][2],float elements[N][8],int i){
    float length = 0;
    float del_x = nodes[(int)elements[i][1]-1][0]-nodes[(int)elements[i][0]-1][0];
    float del_y = nodes[(int)elements[i][1]-1][1]-nodes[(int)elements[i][0]-1][1];
    length = sqrt(pow(del_x,2)+pow(del_y,2));
    elements[i][2]=length;
    return length;
}//End of length

//Compute direction cosines
void cosines(int nodes[N][2],float elements[N][8],int i){
    float L=0,M=0,length=0;
    float del_x = nodes[(int)elements[i][1]-1][0]-nodes[(int)elements[i][0]-1][0];
    float del_y = nodes[(int)elements[i][1]-1][1]-nodes[(int)elements[i][0]-1][1];
    length = Length(nodes,elements,i);
    L = del_x/length;
    M = del_y/length;
    elements[i][6] = L;
    elements[i][7] = M;
}//End of cosines

//Multiplication of Matrix
void matrix_mul(float a[6][6],float b[6][6],float c[6][6],int n,int m){
    int i,j,k;
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            c[i][j]=0;
        }
    }
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
           for(k=0;k<n;k++){
            c[i][j]+=a[i][k]*b[k][j];
           }
        }
    }
}//End of matrix_mul

void mmult(float a[K][2*K],float b[K][1],float c[K][1],int n,int m){
    int i,j,k;
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            c[i][j]=0;
        }
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            c[i][0]+=a[i][j]*b[j][0];
        }
    }
}

void mmult2(float a[K][K],float b[K][1],float c[K][1],int n,int m){
    int i,j,k;
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            c[i][j]=0;
        }
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            c[i][0]+=a[i][j]*b[j][0];
        }
    }
}

//Transpose of Matrix
void transpose(float a[6][6],float b[6][6]){
    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            b[i][j]=a[j][i];
        }
    }
}

//Generate Tranformation matrix
void transformation_matrix(float transformation[6][6],float elements[N][8],int e){
    float l = elements[e][6];
    float m = elements[e][7];
    float transformation_temp[6][6]={{l,m,0,0,0,0},{(-1*m),l,0,0,0,0},{0,0,1,0,0,0},{0,0,0,l,m,0},{0,0,0,(-1*m),l,0},{0,0,0,0,0,1}};
    for(int i=0;i<6;i++)
        for(int j=0;j<6;j++)
            transformation[i][j]=transformation_temp[i][j];
}//End of transformation_matrix

//Assembly of Global Stiffness Matrix
void global_stiffness_assembly(float elements[N][8],int e,int total_elements,int n){
    int p = 3*elements[e][0]-3;
    int q = 3*elements[e][1]-3;
    int i,j;
    for(i=0;i<6;i++){
        for(j=0;j<6;j++){
            if(i<3&&j<3){
                global_stiffness[p+i][p+j]+=element_stiffness[i][j];
            }
            else
            if(i<3&&j>=3){
                global_stiffness[p+i][q+j-3]+=element_stiffness[i][j];
            }
            else
            if(i>=3&&j<3){
                global_stiffness[q+i-3][p+j]+=element_stiffness[i][j];
            }
            else{
                global_stiffness[q+i-3][q+j-3]+=element_stiffness[i][j];
            }
        }
    }
    cout<<"\nGLOBAL STIFFNESS MATRIX AFTER ELEMENT "<<e+1<<endl;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            cout<<setw(10)<<setprecision(3)<<global_stiffness[i][j]<<"|";
        }
        cout<<endl;
    }
}

void print(float a[6][6]){
    int i,j;
    for(i=0;i<6;i++){
        for(j=0;j<6;j++){
            cout<<setw(10)<<setprecision(3)<<a[i][j];
        }
        cout<<endl;
    }
    getch();
}
//Generate Element Stiffness Matrix
void elementstiffnesscalc(float elements[N][8],int e,int total_elements,int n){
    float Le = elements[e][2];
    float E = elements[e][3];
    float A = elements[e][4];
    float I = elements[e][5];
    float element_stiffness_temp[6][6]={{(E*A)/Le,0,0,(-1*E*A)/Le,0,0},
                                  {0,(12*E*I/pow(Le,3)),(6*E*I/pow(Le,2)),0,(-12*E*I/pow(Le,3)),(6*E*I/pow(Le,2))},
                                  {0,(6*E*I/pow(Le,2)),(4*E*I/Le),0,(-6*E*I/pow(Le,2)),(2*E*I/Le)},
                                  {(-1*E*A/Le),0,0,(E*A/Le),0,0},
                                  {0,(-12*E*I/pow(Le,3)),(-6*E*I/pow(Le,2)),0,(12*E*I/pow(Le,3)),(-6*E*I/pow(Le,2))},
                                  {0,(6*E*I/pow(Le,2)),(2*E*I)/Le,0,(-6*E*I/pow(Le,2)),(4*E*I/Le)}
                                 };
    float transformation[6][6];
    float transformation_transpose[6][6];
    transformation_matrix(transformation,elements,e);
    //print(transformation);
    transpose(transformation,transformation_transpose);
    //print(element_stiffness_temp);
    //print(transformation_transpose);
    float temp[6][6];
    matrix_mul(transformation_transpose,element_stiffness_temp,temp,6,6);
    matrix_mul(temp,transformation,element_stiffness,6,6);
    cout<<endl;
    cout<<"ELEMENT STIFFNESS MATRIX"<<e+1<<endl;
    print(element_stiffness);
    global_stiffness_assembly(elements,e,total_elements,3*n);
    for(int i=0;i<K;i++){
        for(int j=0;j<K;j++){
            modified_stiffnesss[i][j]=global_stiffness[i][j];
        }
    }
}//End of elementstiffnesscalc


//Define Elements
void Elements(float elements[N][8],int nodes[N][2],int *e,int n){
    int i = *e;
    int j = 0;
    if(i==0){
        system("cls");
        cout<<"No element exist!"<<endl;
    }
    else{
        View(elements,i,n);
    }
    char ch;
    cout<<"\n\n";
    cout<<"Press B to go back. Enter to Add an Element";
    ch = getch();
    if(ch == 'B'||ch == 'b'){
        return;
    }
    cout<<endl;
    cout<<"Enter the values for element "<<i+1<<endl;
    cout<<"Enter Node i :"; cin>>elements[i][0];
    cout<<"Enter Node j :"; cin>>elements[i][1];
    if(elements[i][1] == elements[i][0]){
        cout<<"Element should be between two distinct nodes.";
    }
    else
    if(checkelement(elements,i,n)==0){
        cout<<"Sorry ! Node doesn't exist.";
    }
    else{
        cosines(nodes,elements,i);
        cout<<"Enter the Young's Modulus of the element (in Pa)"; cin>>elements[i][3];
        cout<<"Enter the Area of the element (in m^2)"; cin>>elements[i][4];
        cout<<"Enter the Moment of Inertia of the element (in m^4)"; cin>>elements[i][5];
        cout<<"Element Successfully Added";
        *e=i+1;
        elementstiffnesscalc(elements,i,i+1,n);
        getch();
    }

    cout<<"Press any key to go back.";
    getch();


}//End of Elements

//Display Global Force Vector
void view_force(int n){
    cout<<"\n   GLOBAL FORCE VECTOR"<<endl;
    for(int i=0;i<3*n;i++){
        cout<<global_force[i][0]<<endl;
    }
}//End of view_force

//Add Forces
void Forces(float elements[N][8],int n,int e){
    system("cls");
    int i,j;
    int type = 0;
    int node;
    int node_i;
    int node_j;
    int element;
    int force_type=0;
    float alpha;
    float delT;
    float magnitude;
    cout<<"Enter the type of loading (1-> Point, 2-> Distributed, 3->Temperature(FOR BARS ALONE) "; cin>>type;
    if(type == 1){
        cout<<"Enter the node at which you want to apply "; cin>>node;
        if(node>n){
            cout<<"Such a node doesn't exist. Enter again";
            getch();
            return;
        }
        cout<<"Enter the whether force/moment (1->Rz 2->Fy 3->Fx) "; cin>>force_type;
        if(force_type >3 || force_type <1 ){
            cout<<"Select a valid option. (1-3)";
            getch();
            return;
        }
        cout<<"Enter the magnitude(N) "; cin>>magnitude;
        global_force[(3*node)-(force_type)][0]+=magnitude;
    }
    else
    if(type == 2){
        cout<<"Enter the element on which it is acting "; cin>>element;
        if(element>e){
            cout<<"Such a element doesn't exist. Enter again";
            getch();
            return;
        }
        cout<<"Enter the magnitude of distributed force(N/m)"; cin>>magnitude;
        node_i = elements[element-1][0];
        node_j = elements[element-1][1];
        global_force[(3*node_i)-3][0]+=(elements[element-1][7]*magnitude*elements[element-1][2])/2;
        global_force[(3*node_i)-2][0]+=(elements[element-1][6]*magnitude*elements[element-1][2])/2;
        global_force[(3*node_i)-1][0]+=(magnitude*pow(elements[element-1][2],2))/12;
        global_force[(3*node_j)-3][0]+=(elements[element-1][7]*magnitude*elements[element-1][2])/2;
        global_force[(3*node_j)-2][0]+=(elements[element-1][6]*magnitude*elements[element-1][2])/2;
        global_force[(3*node_j)-1][0]+=(-1*magnitude*pow(elements[element-1][2],2))/12;
    }
    else
    if(type == 3){
        cout<<"Enter the element subject to temperature effect "; cin>>element;
        if(element>e){
            cout<<"Such a element doesn't exist. Enter again";
            getch();
            return;
        }
        node_i = elements[element-1][0];
        node_j = elements[element-1][1];

        cout<<"Enter the Coeffecient of Linear Thermal Expansion of the element (/C)"; cin>>alpha;
        cout<<"Enter the temperature gradient (C) "; cin>>delT;
        global_force[3*node_i-3][0]+=(elements[element-1][3]*elements[element-1][4]*alpha*delT*(-1*elements[element-1][6]));
        global_force[3*node_i-2][0]+=(elements[element-1][3]*elements[element-1][4]*alpha*delT*(-1*elements[element-1][7]));
        global_force[3*node_j-3][0]+=(elements[element-1][3]*elements[element-1][4]*alpha*delT*(elements[element-1][6]));
        global_force[3*node_j-2][0]+=(elements[element-1][3]*elements[element-1][4]*alpha*delT*(elements[element-1][7]));
    }
    else{
        cout<<"Select a valid option. (1-3)";
        getch();
        return;
    }
    cout<<"Force Successfully added. Enter to view Force vector";
    for(j=0; j<3*n;j++){
        modified_force[j][0]=global_force[j][0];
    }
    getch();
    view_force(n);
    getch();
}//End of Forces

//Display Global Displacement Vector
void view_displacement(int n){
    cout<<"\nGLOBAL DISPLACEMENT VETOR "<<endl;
    for(int i=0;i<3*n;i++){
        if(global_displacement[i][1]){
           cout<<global_displacement[i][0]<<endl;
        }
        else{
            cout<<"q"<<i+1  <<endl;
        }
    }
}//End of view_displacement

//Define Boundary Conditions
void Boundary_conditions(float elements[N][8],int n,int e){
    system("cls");
    cout<<"BOUNDARY CONDITIONS"<<endl;
    float C = 0;
    int i,j,k;
    for(i=0;i<3*n;i++){
        for(j=0;j<3*n;j++){
            if(global_stiffness[i][j]>C)
                C = global_stiffness[i][j];
        }
    }
    C *= 10000;
    int type;
    int node;
    float value = 0;
    int dof;
    int a,b,c,Qp1,Qp2;
    cout<<"Select the type of constraint. (1-> Single Point. 2-> Multi Point) "; cin>>type;
    if(type == 1){
        cout<<"Enter the DoF which you know (1 - "<<(3*n)<<") "; cin>>dof;
        if(dof>3*n||dof<1){
            cout<<"Such a node doesn't exist. Enter again";
            getch();
            return;
        }
        if(dof%3==0)
            cout<<"Enter known rotation ";
        else
            cout<<"Enter known displacement ";
        cin>>value;
        global_displacement[dof-1][1]=1;
        global_displacement[dof-1][0]=value;
        modified_stiffnesss[dof-1][dof-1]+=C;
        modified_force[dof-1][0]+=(C*value);
    }
    else
    if(type == 2){
        cout<<"MPC is of the form : a*Qp1 + b*Qp2 = c"<<endl;
        cout<<"Enter a "; cin>>a;
        cout<<"Enter Qp1 ";cin>>Qp1;
        cout<<"Enter b ";cin>>b;
        cout<<"Enter Qp2 ";cin>>Qp2;
        cout<<"Enter c ";cin>>c;
        for(i=0;i<3*n;i++){
            for(j=0;j<3*n;j++){
                if(i==j&&i==Qp1){
                    modified_stiffnesss[i][j]+=(C*pow(a,2));
                    modified_force[i][0]+=(C*c*a);
                }
                else
                if((i==Qp1 && j==Qp2)||(i==Qp2 && j == Qp1)){
                    modified_stiffnesss[i][j]+=(C*a*b);
                }
                else
                if(i==j&&i==Qp2){
                    modified_stiffnesss[i][j]+=(C*pow(b,2));
                    modified_force[i][0]+=(C*c*b);
                }
            }
        }
    }
    cout<<"Constraint Successfully applied. Enter to view the displacement vector";
    view_displacement(n);
    getch();
}//End of Boundary_conditions

//Returns absolute value
float mod(float n){
    if(n<0){
        return (-1*n);
    }
    else{
        return n;
    }
}//End of mod

//View K_Mod
void view_modified_stiffness(int n){
    for(int i=0;i<3*n;i++){
        for(int j=0;j<3*n;j++){
            cout<<setw(10)<<setprecision(3)<<modified_stiffnesss[i][j];
        }
        cout<<endl;
    }
}//End of view_modified_stiffness

//View F_Mod
void view_modified_force(int n){
    for(int i=0;i<3*n;i++){
        cout<<modified_force[i][0]<<endl;
    }
}//End of view_modified_force

//Compute the inverse of stiffness matrix
void inverse(float a[K][2*K],int n){

   int i,j,k;
   float t;
   for(i=0;i<n;i++)
   {
      for(j=n;j<2*n;j++)
      {
          if(i==j-n)
             a[i][j]=1;
         else
             a[i][j]=0;
       }
   }
   for(i=0;i<n;i++)
   {
      t=a[i][i];
      for(j=i;j<2*n;j++)
          a[i][j]=a[i][j]/t;
      for(j=0;j<n;j++)
      {
         if(i!=j)
         {
            t=a[j][i];
            for(k=0;k<2*n;k++)
                a[j][k]=a[j][k]-t*a[i][k];
          }
      }
   }
   system("cls");
   for(i=0;i<n;i++){
        for(j=n;j<2*n;j++){
            a[i][j-n]=a[i][j];
            cout<<setw(10
                        )<<setprecision(3)<<a[i][j]<<"|";
        }
        cout<<endl;
   }
   getch();
}//End of inverse

//Solve using FEM
void Solve(float elements[N][8],int e,int n){
    system("cls");
    int i,j,k;
    cout<<"MODIFIED STIFFNESS MATRIX"<<endl;
    view_modified_stiffness(n);
    cout<<"MODIFIED FORCE VECTOR"<<endl;
    view_modified_force(n);
    getch();
    float modified_stiffnesss_inverse[K][2*K];
    for(i=0;i<3*n;i++){
        for(j=0;j<3*n;j++){
            modified_stiffnesss_inverse[i][j]=modified_stiffnesss[i][j];
        }
    }
    inverse(modified_stiffnesss_inverse,3*n);

    mmult(modified_stiffnesss_inverse,modified_force,solution,3*n,1);

    system("cls");
    cout<<"NODAL SOLUTION"<<endl;
    for(i=0;i<3*n;i++){
        cout<<"q "<<i+1<<" : "<<solution[i][0]<<endl;
    }

    getch();
}//end of Solve

//Calculation and Display of Secondary Variable
void SecondaryVariable(float elements[N][8],int e,int n){
    float kq[K][1];
    mmult2(global_stiffness,solution,kq,3*n,1);
    float reaction_force[K][1];
    int i;
    for(i=0;i<3*n;i++){
        reaction_force[i][0] = kq[i][0]-global_force[i][0];
    }
    system("cls");
    cout<<"REACTION FORCES "<<endl;
    for(i=0;i<3*n;i++){
        cout<<reaction_force[i][0]<<endl ;
    }
    getch();
}
//End of SecondaryVariable

//Thank you note
void Quit(){
    system("cls");
    cout<<"THANK YOU FOR USING THIS SOLVER";
    cout<<endl<<"PROJECT BY ";
    cout<<endl<<"S.Sidarth";
    cout<<endl<<"MDM16B032";
}//End of Quit


int main(){

    char ch;
    int n=0,e=0;
    int nodes[N][2];
    float elements[N][8];
    settozero(element_stiffness,global_stiffness,global_force,2);
    settozero(element_stiffness,global_stiffness,global_force,3);
    do{
        ch = menu();
        switch(ch){
            case '1' : { Nodes(nodes,&n); break; }
            case '2' : { Elements(elements,nodes,&e,n); break; }
            case '3' : { View(elements,e,n); getch(); break; }
            case '4' : { Forces(elements,n,e); break; }
            case '5' : { Boundary_conditions(elements,n,e); break; }
            case '6' : { Solve(elements,e,n); break; }
            case '7' : { SecondaryVariable(elements,e,n); break; }
            case '8' : { Quit(); return(0); }
            default : { cout<<"\nWrong Choice. Please Enter between 1-6. Press any key to go enter again."; getch(); break; }
        }
    }while(ch!='8');

    return 0;
}//End of Program
