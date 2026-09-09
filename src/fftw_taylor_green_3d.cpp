#include <fftw3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Complex = std::complex<double>;
using RealField = std::vector<double>;
using SpectralField = std::vector<Complex>;
using RealVectorField = std::array<RealField, 3>;
using SpectralVectorField = std::array<SpectralField, 3>;

struct Grid {
    int n;
    std::size_t size() const { return static_cast<std::size_t>(n) * n * n; }
    std::size_t index(int i, int j, int k) const { return (static_cast<std::size_t>(i) * n + j) * n + k; }
    int wave_number(int q) const { return (q <= n / 2) ? q : q - n; }
};

class Fft3d {
public:
    explicit Fft3d(const Grid& grid) : grid_(grid) {
        input_ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * grid.size()));
        output_ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * grid.size()));
        if (!input_ || !output_) throw std::runtime_error("FFTW allocation failed");
        forward_ = fftw_plan_dft_3d(grid.n, grid.n, grid.n, input_, output_, FFTW_FORWARD, FFTW_ESTIMATE);
        backward_ = fftw_plan_dft_3d(grid.n, grid.n, grid.n, input_, output_, FFTW_BACKWARD, FFTW_ESTIMATE);
        if (!forward_ || !backward_) throw std::runtime_error("FFTW plan failed");
    }
    ~Fft3d() {
        if (forward_) fftw_destroy_plan(forward_);
        if (backward_) fftw_destroy_plan(backward_);
        fftw_free(input_); fftw_free(output_);
    }
    SpectralField forward(const RealField& values) {
        for (std::size_t i=0;i<grid_.size();++i){input_[i][0]=values[i];input_[i][1]=0.0;}
        fftw_execute(forward_);
        const double s=1.0/static_cast<double>(grid_.size());
        SpectralField out(grid_.size());
        for(std::size_t i=0;i<grid_.size();++i) out[i]=s*Complex{output_[i][0],output_[i][1]};
        return out;
    }
    RealField inverse(const SpectralField& spec) {
        for(std::size_t i=0;i<grid_.size();++i){input_[i][0]=spec[i].real();input_[i][1]=spec[i].imag();}
        fftw_execute(backward_);
        RealField out(grid_.size());
        for(std::size_t i=0;i<grid_.size();++i) out[i]=output_[i][0];
        return out;
    }
private:
    Grid grid_;
    fftw_complex* input_{nullptr};
    fftw_complex* output_{nullptr};
    fftw_plan forward_{nullptr};
    fftw_plan backward_{nullptr};
};

SpectralVectorField make_spectral(std::size_t n){ SpectralVectorField f; for(auto& c:f)c.assign(n,{0,0}); return f; }
RealVectorField make_real(std::size_t n){ RealVectorField f; for(auto& c:f)c.assign(n,0.0); return f; }

void project_dealias(const Grid& g,SpectralVectorField& f){
    const int cutoff=g.n/3;
    for(int i=0;i<g.n;++i){int kx=g.wave_number(i);for(int j=0;j<g.n;++j){int ky=g.wave_number(j);for(int k=0;k<g.n;++k){int kz=g.wave_number(k);auto idx=g.index(i,j,k);
        if(std::abs(kx)>cutoff||std::abs(ky)>cutoff||std::abs(kz)>cutoff){for(auto& c:f)c[idx]={0,0};continue;}
        double k2=double(kx*kx+ky*ky+kz*kz); if(k2==0)continue;
        Complex kd=double(kx)*f[0][idx]+double(ky)*f[1][idx]+double(kz)*f[2][idx];
        f[0][idx]-=double(kx)*kd/k2; f[1][idx]-=double(ky)*kd/k2; f[2][idx]-=double(kz)*kd/k2;
    }}}
}

SpectralVectorField omega_hat(const Grid& g,const SpectralVectorField& u){
    SpectralVectorField w=make_spectral(g.size()); const Complex I{0,1};
    for(int i=0;i<g.n;++i){double kx=g.wave_number(i);for(int j=0;j<g.n;++j){double ky=g.wave_number(j);for(int k=0;k<g.n;++k){double kz=g.wave_number(k);auto idx=g.index(i,j,k);
        w[0][idx]=I*(ky*u[2][idx]-kz*u[1][idx]);
        w[1][idx]=I*(kz*u[0][idx]-kx*u[2][idx]);
        w[2][idx]=I*(kx*u[1][idx]-ky*u[0][idx]);
    }}}
    return w;
}

SpectralVectorField rhs(const Grid& g,Fft3d& fft,const SpectralVectorField& u,double nu){
    auto wh=omega_hat(g,u); RealVectorField ur=make_real(g.size()),wr=make_real(g.size());
    for(int c=0;c<3;++c){ur[c]=fft.inverse(u[c]);wr[c]=fft.inverse(wh[c]);}
    RealVectorField rot=make_real(g.size());
    for(std::size_t q=0;q<g.size();++q){rot[0][q]=ur[1][q]*wr[2][q]-ur[2][q]*wr[1][q];rot[1][q]=ur[2][q]*wr[0][q]-ur[0][q]*wr[2][q];rot[2][q]=ur[0][q]*wr[1][q]-ur[1][q]*wr[0][q];}
    SpectralVectorField out=make_spectral(g.size()); for(int c=0;c<3;++c)out[c]=fft.forward(rot[c]); project_dealias(g,out);
    for(int i=0;i<g.n;++i){int kx=g.wave_number(i);for(int j=0;j<g.n;++j){int ky=g.wave_number(j);for(int k=0;k<g.n;++k){int kz=g.wave_number(k);auto idx=g.index(i,j,k);double k2=double(kx*kx+ky*ky+kz*kz);for(int c=0;c<3;++c)out[c][idx]-=nu*k2*u[c][idx];}}}
    return out;
}

SpectralVectorField add_scaled(const SpectralVectorField&a,const SpectralVectorField&b,double s){auto o=a;for(int c=0;c<3;++c)for(std::size_t i=0;i<o[c].size();++i)o[c][i]+=s*b[c][i];return o;}
SpectralVectorField rk4(const Grid&g,Fft3d&fft,const SpectralVectorField&u,double dt,double nu){auto k1=rhs(g,fft,u,nu);auto k2=rhs(g,fft,add_scaled(u,k1,.5*dt),nu);auto k3=rhs(g,fft,add_scaled(u,k2,.5*dt),nu);auto k4=rhs(g,fft,add_scaled(u,k3,dt),nu);auto o=u;for(int c=0;c<3;++c)for(std::size_t i=0;i<o[c].size();++i)o[c][i]+=(dt/6.0)*(k1[c][i]+2.0*k2[c][i]+2.0*k3[c][i]+k4[c][i]);project_dealias(g,o);return o;}

RealVectorField initial_3d_tg(const Grid& g){
    RealVectorField u=make_real(g.size()); const double tp=2.0*std::numbers::pi;
    for(int i=0;i<g.n;++i){double x=tp*i/g.n;for(int j=0;j<g.n;++j){double y=tp*j/g.n;for(int k=0;k<g.n;++k){double z=tp*k/g.n;auto idx=g.index(i,j,k);u[0][idx]=std::sin(x)*std::cos(y)*std::cos(z);u[1][idx]=-std::cos(x)*std::sin(y)*std::cos(z);u[2][idx]=0.0;}}}
    return u;
}

void write_snapshot(std::ofstream& out,const Grid& g,Fft3d& fft,const SpectralVectorField& uh,double t){
    RealVectorField u=make_real(g.size()); auto wh=omega_hat(g,uh); RealVectorField w=make_real(g.size());
    for(int c=0;c<3;++c){u[c]=fft.inverse(uh[c]);w[c]=fft.inverse(wh[c]);}
    const double tp=2.0*std::numbers::pi; out<<std::setprecision(17);
    for(int plane=0;plane<3;++plane){int fixed=g.n/2;for(int a=0;a<g.n;++a){for(int b=0;b<g.n;++b){int i=0,j=0,k=0;if(plane==0){i=a;j=b;k=fixed;}else if(plane==1){i=a;j=fixed;k=b;}else{i=fixed;j=a;k=b;}auto idx=g.index(i,j,k);double wm=std::sqrt(w[0][idx]*w[0][idx]+w[1][idx]*w[1][idx]+w[2][idx]*w[2][idx]);out<<t<<','<<plane<<','<<tp*i/g.n<<','<<tp*j/g.n<<','<<tp*k/g.n<<','<<u[0][idx]<<','<<u[1][idx]<<','<<u[2][idx]<<','<<w[0][idx]<<','<<w[1][idx]<<','<<w[2][idx]<<','<<wm<<'\n';}}}
}

}

int main(int argc,char**argv){
    int n=24; double final_time=0.5,dt=0.002,nu=0.01; std::string out_path="taylor_green_3d_snapshots.csv"; int snapshots=21;
    if(argc>1)n=std::stoi(argv[1]); if(argc>2)final_time=std::stod(argv[2]); if(argc>3)dt=std::stod(argv[3]); if(argc>4)nu=std::stod(argv[4]); if(argc>5)out_path=argv[5]; if(argc>6)snapshots=std::stoi(argv[6]);
    if(n<8||n%2!=0||final_time<=0||dt<=0||nu<=0||snapshots<2)throw std::invalid_argument("invalid arguments");
    int steps=static_cast<int>(std::ceil(final_time/dt)); dt=final_time/steps; Grid g{n}; Fft3d fft(g); auto initial=initial_3d_tg(g); auto uh=make_spectral(g.size()); for(int c=0;c<3;++c)uh[c]=fft.forward(initial[c]); project_dealias(g,uh);
    std::ofstream out(out_path); if(!out)throw std::runtime_error("failed to open output"); out<<"time,plane,x,y,z,u,v,w,omega_x,omega_y,omega_z,omega_mag\n";
    int next_snap=0;
    for(int step=0;step<=steps;++step){
        double t=step*dt;
        while(next_snap<snapshots){
            int target_step=static_cast<int>(std::llround(static_cast<double>(next_snap)*steps/(snapshots-1)));
            if(step<target_step)break;
            write_snapshot(out,g,fft,uh,t);
            ++next_snap;
        }
        if(step<steps)uh=rk4(g,fft,uh,dt,nu);
    }
    std::cout<<"wrote 3D Taylor-Green snapshots: "<<out_path<<"\n"; return 0;
}
