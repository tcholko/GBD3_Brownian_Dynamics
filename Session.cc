#include "Session.h"
#include "Model.h"
#include "Body.h"
#include "BindingCriteria.h"
#include "Grid_EX.h"
#include "Grid.h"
#include "Strings.h"
#include "Grid_ES.h"
#include "Grid_D.h"
#include "Timer.h"


///Session definition 
Session::Session(Model *m, SimulationConfig s) {
  id = 0;
  model = m;
  type = s;
  Nreplicates = 0;
  Nbind.set_value(0);
  Nexit.set_value(0);
  Ntlim.set_value(0);
  Davg = 0.;
  done = false;
}

Session::~Session() {
}


void Session::populateLigands() {
  for(int i=0; i < Nreplicates; i++) {
    int rndConf = floor(random(0.0, (double)conformations.size())); //Randomly select a ligand conformation
    Body *rc = conformations[rndConf];
    cout << "about to create a new Body\n for replicate " << i+1<< endl;
    Body *bi = new Body(model, this); // Construct new Body object by passing parameters "model" and "this"
    for(int j=0; j < rc->beads.size(); j++) { //Copy over every bead in conformations to bi
      Bead *nb = new Bead();
      Bead *rb = rc->beads[j]; //*rb is pointer to jth bead in the conformation
      nb->R.x = rb->R.x;       //Next, coords and parameters of *rb are copied to *nb
      nb->R.y = rb->R.y;
      nb->R.z = rb->R.z;
      nb->q = rb->q;
      nb->r = rb->r;
      nb->m = rb->m;
      nb->type = rc->beads[j]->type;
      bi->beads.push_back(nb); //*bi is the Body (i.e. ligand) made of all these beads (i.e. atoms). Each Bead nb is pushed into
    }                          //Body bi's bead vector
    bi->define();
    if(model->useRestart) { 
      positionLigandRestart(bi, i); 
      cout << "using rst time from vrt[" << i << "]\n";
      bi->t = model->vrt[i];
      cout<< "t for ligand "<<i+1<< " is " << model->vrt[i] <<endl;
     }
    else { positionLigand(bi); }   
    model->ligands.push_back(bi);   //push ligand bi into vector holding all ligands accross all sessions
    ligands.push_back(bi);          //push ligand into vector holding ligands in one session instance
  }

  Davg = 0.;
  for(int i=0; i < conformations.size(); i++) {
    Davg += conformations[i]->D;
  }
  Davg /= conformations.size();
}


void Session::recordBeta() {
  double _Nbind = (double)Nbind.get_value();
  double _Ndone = (double)Nexit.get_value() + (double)Ntlim.get_value();
  double beta = _Nbind / (_Nbind + _Ndone);
  beta_history.push_back(beta);
  while(beta_history.size() > model->convergence_window) {
    beta_history.pop_front();
  }
}


double Session::checkConvergence() {
  if(beta_history.size() < model->convergence_window) return -1;

  double m = 0., s = 0.;
  for(int i=0; i < beta_history.size(); i++) {
    m += beta_history[i];
  }
  m /= beta_history.size();
  if(m == 0) return -1;

  for(int i=0; i < beta_history.size(); i++) {
    s += pow(beta_history[i]-m, 2);
  }
  s = sqrt(s / beta_history.size());
  double soverm = s / m;
  if(soverm <= model->convergence and soverm != 0. and done == false) {
    model->lout << "* Convergence criteria reached. Exiting successfully." << endl;
    done = true;
    for(int i=0; i < ligands.size(); i++) ligands[i]->done = true;
  }

  return soverm;
}




//SessionNAM definition
SessionNAM::SessionNAM(Model *m) : Session(m, CONFIGURATION_RADIAL) {
  b = 0.;
  q = 0.;
  q2 = 0.;
}

void SessionNAM::positionLigandRestart(Body *bi, int iter) {
 vertex Q;
 int bcnt=0;
 double tqx=0., tqy=0., tqz=0., tm=0.;
  bi->center(); // Center before changing R.x, etc., center() applies an equal/opp vector to Body's COM
  for (int r=bi->beads.size()*3*iter; r<bi->beads.size()*3*(iter+1); r+=3){
     bi->beads[bcnt]->R.x = model->vrc[r];
     bi->beads[bcnt]->R.y = model->vrc[r+1];
     bi->beads[bcnt]->R.z = model->vrc[r+2];
     tqx += model->vrc[r]  *  bi->beads[bcnt]->m ;
     tqy += model->vrc[r+1] * bi->beads[bcnt]->m ;
     tqz += model->vrc[r+2] * bi->beads[bcnt]->m ;
     tm += bi->beads[bcnt]->m;
     //cout << "bead " << bcnt<< " coords = "<<bi->beads[bcnt]->R.x<<" "<<bi->beads[bcnt]->R.y<<" "<<bi->beads[bcnt]->R.z<<endl;
     bcnt++;
    }
  Q.x = tqx/tm; Q.y=tqy/tm; Q.z=tqz/tm;
  cout<<"ligdn COM is " <<Q.x<<" "<<Q.y<<" "<<Q.z<<endl;
  bi->R.x = Q.x; //Body COM must be moved too or the dynamics are super fast and disobey boundaries. Still not sure why
  bi->R.y = Q.y; //but I know it's caused by moving beads but not moving the COM with them
  bi->R.z = Q.z;
}

void SessionNAM::positionLigand(Body *bi) {
 if (model->rand_start == 0 && b > model->exBnd){ 
  cout << "! Error: Ligand starting positions are outside of simulation boundary! (must have b <= exbound)" << endl; 
  std::exit(0); 
  }
 cout << "inside sessionNAM postionLigand\nstep = " << model->step <<endl;

vertex Q = { random(-1., 1.), random(-1., 1.), random(-1., 1.) };
// Create starting square plane
 if (model->boundary == 1) {
  Q.x *= b;
  Q.y *= b;
  Q.z = model->bPlnHgt;
 
   if (model->rand_start == 1) {
     bool penetrating = true;
      while (penetrating) {
       Q.x  = random(-1*model->exBnd, model->exBnd);
       Q.y  = random(-1*model->exBnd, model->exBnd);
       Q.z  = random(0., 500.); // Change later !!!
       int pen_cnt = 0;
        for (int i=0; i<model->x_pen_chk.size(); i++){
          double x_pen2 = pow((model->x_pen_chk[i] - Q.x), 2);
          double y_pen2 = pow((model->y_pen_chk[i] - Q.y), 2);
          double z_pen2 = pow((model->z_pen_chk[i] - Q.z), 2);
          double d_pen2 = x_pen2+y_pen2+z_pen2;
           if (d_pen2 > 100.) { pen_cnt++; } 
           if (pen_cnt >= model->x_pen_chk.size()) { penetrating = false; }
        }
      }
    }
  }
// Or create starting hexagonal plane
 else if (model->boundary == 2) {
  Q.y *= b; // Q.y is some % of b (input file), then Q.x is a little larger for slanted sides of hexagon 
  if (Q.y >= 0) {Q.x *= (b-Q.y)*tan(M_PI/6) + b*tan(M_PI/6);}
      else      {Q.x *= (-b-Q.y)*tan(M_PI/6) - b*tan(M_PI/6);}
  Q.z = model->bPlnHgt;

    if (model->rand_start == 1) {
      bool penetrating = true;
       while (penetrating) {
        Q.y = random(-1*model->exBnd, model->exBnd);
         if (Q.y >= 0) { Q.x = random(-1 * ((model->exBnd-Q.y)*tan(M_PI/6) + model->exBnd*tan(M_PI/6)), ((model->exBnd-Q.y)*tan(M_PI/6) + model->exBnd*tan(M_PI/6))); }
          else        { Q.x = random(-1 * ((-model->exBnd-Q.y)*tan(M_PI/6) - model->exBnd*tan(M_PI/6)), ((-model->exBnd-Q.y)*tan(M_PI/6) - model->exBnd*tan(M_PI/6))); }
        Q.z = random(0., 500.); // Make this between floor and ceiling variables later
        int pen_cnt=0;
         for (int i=0; i<model->x_pen_chk.size(); i++){
           double x_pen2 = pow((model->x_pen_chk[i] - Q.x), 2);
           double y_pen2 = pow((model->y_pen_chk[i] - Q.y), 2);
           double z_pen2 = pow((model->z_pen_chk[i] - Q.z), 2);
           double d_pen2 = x_pen2+y_pen2+z_pen2;
            if (d_pen2 > 100.) { pen_cnt++; }
            if (pen_cnt >= model->x_pen_chk.size()) { penetrating = false; }
        }
      }
    }
  }
// Or create spherical plane
 else if (model->boundary == 3) {
  double l = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
  Q.x *= b / l;
  Q.y *= b / l;
  Q.z *= b / l;
  
   if (model->rand_start == 1) {
    bool penetrating = true;
    while (penetrating) {
      bool out_of_bounds=true;
      while(out_of_bounds){
        Q.x = random(-q, q); 
        Q.y = random(-q, q);
        Q.z = random(-q, q);
        double dx2 = pow((Q.x - model->center.x), 2);
        double dy2 = pow((Q.y - model->center.y), 2);
        double dz2 = pow((Q.z - model->center.z), 2);
        if (dx2+dy2+dz2 < pow(q, 2)) { out_of_bounds=false; }
        }
      int pen_cnt = 0;
        for (int i=0; i<model->x_pen_chk.size(); i++){
          double x_pen2 = pow((model->x_pen_chk[i] - Q.x), 2);
          double y_pen2 = pow((model->y_pen_chk[i] - Q.y), 2);
          double z_pen2 = pow((model->z_pen_chk[i] - Q.z), 2);
          double d_pen2 = x_pen2+y_pen2+z_pen2;
          if (d_pen2 > 100.) { pen_cnt++; }
          if (pen_cnt >= model->x_pen_chk.size()) { penetrating = false; }
        }
      }
    }
  }
  bi->center();
  bi->translate(model->center.x + Q.x, model->center.y + Q.y, /*model->center.z +*/ Q.z); //Move ligand into place w/ translate function
  bi->rotate(random(0., M_PI), random(0., M_PI), random(0., M_PI)); //Give it random rotational orientation
  cout << "just positioned a ligand in sessionNAM::positionLigand\n";
}

void SessionNAM::printRateConstant() {
  int Ndone = Nbind.get_value() + Nexit.get_value();
  double conv = checkConvergence(); 

  if(Ndone == 0) {
    double kb = 4. * M_PI * b * Davg   *Na*1e12*1e-27;
    model->lout << "   (session " << id << ") ";
    model->lout << "kd(b)=" << kb << " ";
    model->lout << "b=" << b << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
    return;
  }

  for(int bsi=0; bsi < bindingCriteria.size(); bsi++) {
    double tavg = bindingCriteria[bsi]->t_avgt.get_value() / bindingCriteria[bsi]->Nbind.get_value();
    double B = ((double)bindingCriteria[bsi]->Nbind.get_value()) / ((double)Ndone);
    double kb = 4. * M_PI * b * Davg;
    //double vol = 4*exBnd*exBnd*((bounds_max.z + ceiling) - (bounds_min.z + 40));   //40 is there because of padding. Make this use a padding variable instead of        just 40 later TC 071119 
    //double k = (kb * B) / (1 - ((1 - B)*b/q));
    double k = 50000000. / tavg; // New k with no exit criteria, beta not needed. 50,000,000 A^3 approx vol. Make this calculable later
    k *= Na * 1e12 * 1e-27;
    model->lout << "   (session " << id << " bs " << bsi << ") ";
    model->lout << "k_on = " << k << " M⁻¹s⁻¹ ";
    model->lout << "Nbind=" << bindingCriteria[bsi]->Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    model->lout << "β=" << B << " ";
    model->lout << "conv=" << conv << " ";
    if(done)
      model->lout << "(convergence reached) ";
    model->lout << "kd(b)=" << kb * Na * 1e12 * 1e-27 << " ";
    model->lout << "b=" << b << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << "Boundary type = " << model->boundary; //New
    model->lout << endl;
  }

  double B = ((double)Nbind.get_value()) / ((double)Ndone);
  if(bindingCriteria.size() > 1) {
    double tavg = t_avgt.get_value() / Nbind.get_value();
    double kb = 4. * M_PI * b * Davg;
    //double vol = 4*exBnd*exBnd*((bounds_max.z + ceiling) - (bounds_min.z + 40));   //40 is there because of padding. Make this use a padding variable instead of jus    t 40 later TC 071119  
    //double k = (kb * B) / (1 - ((1 - B)*b/q));  
    double k = 50000000. / tavg;    //New k with no exit criteria, beta not needed TC0710
    k *= Na * 1e12 * 1e-27;
    model->lout << "   (session " << id << ") ";
    model->lout << "k_on = " << k << " M⁻¹s⁻¹ ";  
    model->lout << "Nbind=" << Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    model->lout << "β=" << B << " ";
    model->lout << "kd(b)=" << kb * Na * 1e12 * 1e-27 << " ";
    model->lout << "b=" << b << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
  }

  if(model->max_simulations > 0 and Ndone >= model->max_simulations) {
    model->lout << "> Maximum simulations reached. Exiting." << endl;
    model->done = true;
  }
}


void SessionNAM::checkLigand(Body *bi) {
// Only check q sphere penetration for spherical simulation cell
 if (model->boundary == 3){
  double dr[3], l2;
  dr[0] = bi->R.x - model->center.x;
  dr[1] = bi->R.y - model->center.y;
  dr[2] = bi->R.z - model->center.z;
  l2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
  //If distance from center is greater than/equal to q radius then set as exit event and place a new ligand using the
  //same pointer
  if(l2 >= q2) {
    if(model->logExiters)
      model->lout << "#" << id << "\t Escape event at t=" << bi->t << " ps  (t_dwell=" << bi->t_dwell << "ps, max=" << bi->t_dwell_max << "ps, total=" << bi->t_dwell_total << "ps)" << endl;
    bi->exited = true;
    bi->session->positionLigand(bi);
    bi->t = 0.;
    bi->t_dwell = 0.;
    bi->t_dwell_max = 0.;
    bi->t_dwell_total = 0.;
  }
 }
}


//SessionDirect definition
SessionDirect::SessionDirect(Model *m) : Session(m, CONFIGURATION_ABSOLUTE_RADIAL) {
  t_avgt.set_value(0.);
}


void SessionDirect::positionLigand(Body *bi) {
  bi->center();
  bi->translate(start.x, start.y, start.z);
}


void SessionDirect::printRateConstant() {
  int Ndone = Nbind.get_value() + Nexit.get_value();
  if(Ndone == 0) return;
  double conv = checkConvergence();

  for(int bsi=0; bsi < bindingCriteria.size(); bsi++) {
    double B = ((double)bindingCriteria[bsi]->Nbind.get_value()) / ((double)Ndone);
    double tavg = bindingCriteria[bsi]->t_avgt.get_value() / bindingCriteria[bsi]->Nbind.get_value();
    double k = B / (tavg * 1e-12);
    //double k = 1 / (tavg * 1e-12); // New k with no exit criteria, beta not needed TC0710
    
    model->lout << "   (session " << id << " bs " << bsi << ") ";
    model->lout << "k_direct = " << k << " s⁻¹ ";
    model->lout << "Nbind=" << bindingCriteria[bsi]->Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    model->lout << "t_avg=" << tavg << " ";
    model->lout << "β=" << B << " ";
    model->lout << "conv=" << conv << " ";
    if(done)
      model->lout << "(convergence reached) ";
    model->lout << "b=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
  }

  double B = ((double)Nbind.get_value()) / ((double)Ndone);
  if(bindingCriteria.size() > 1) {
    double tavg = t_avgt.get_value() / Nbind.get_value();
    double k = B / (tavg * 1e-12);
    //double k = 1 / (tavg * 1e-12); // New k with no exit criteria, beta not needed TC0710

    model->lout << "   (session " << id << ") ";
    model->lout << "k_direct = " << k << " s⁻¹ ";
    model->lout << "Nbind=" << Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    model->lout << "t_avg=" << tavg << " ";
    model->lout << "β=" << B << " ";
    model->lout << "conv=" << conv << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
  }

  if(model->max_simulations > 0 and Ndone >= model->max_simulations) {
    model->lout << "> Maximum simulations reached. Exiting." << endl;
    model->done = true;
  }
}


void SessionDirect::checkLigand(Body *bi) {
  double dr[3], l2;
  dr[0] = bi->R.x - model->center.x;
  dr[1] = bi->R.y - model->center.y;
  dr[2] = bi->R.z - model->center.z;
  l2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

  if(l2 >= q2) {
    //bi->done = true;
    //*Nexit += 1;
    // Record Beta value after exit event
    //bi->session->recordBeta();
    // Reposition ligand
    bi->session->positionLigand(bi);
    bi->exited = true;
    if(model->logExiters)
      model->lout << "#" << id << "\t Escape event at t=" << bi->t << " ps  l=" << sqrt(l2) << " (t_dwell=" << bi->t_dwell << "ps, max=" << bi->t_dwell_max << "ps, total=" << bi->t_dwell_total << "ps)" << endl;
    bi->t = 0.;
    bi->t_dwell = 0.;
    bi->t_dwell_max = 0.;
    bi->t_dwell_total = 0.;
  }
}



