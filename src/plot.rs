#![allow(unused)]

use plotters::{prelude::*, element::PointCollection, coord::ranged3d::{ProjectionMatrixBuilder, ProjectionMatrix}};

type T = f32;

const CLAMP: T = 200.0;
const PLOT_RES: (u32, u32) = (1024, 760);
const PLOT_CAPTION_FONT: (&str, u32) = ("sans", 20);
const PLOT_MARGIN: u32 = 5;
const PLOT_LABEL_AREA_SIZE: u32 = 30;

fn isometric(mut pb: ProjectionMatrixBuilder) -> ProjectionMatrix
{
    pb.yaw = core::f64::consts::FRAC_PI_4;
    pb.pitch = core::f64::consts::FRAC_PI_4;
    pb.scale = 0.7;
    pb.into_matrix()
}

pub fn plot_curves<const N: usize, const M: usize>(
    plot_title: &str, plot_path: &str,
    x: [[T; N]; M],
    mut y: [[T; N]; M]
) -> Result<(), Box<dyn std::error::Error>>
{
    let x_min = x.into_iter().flatten().filter(|x| x.is_finite()).reduce(T::min).unwrap();
    let x_max = x.into_iter().flatten().filter(|x| x.is_finite()).reduce(T::max).unwrap();
    
    let y_min = y.into_iter().flatten().filter(|x| x.is_finite()).reduce(T::min).unwrap();
    let y_max = y.into_iter().flatten().filter(|x| x.is_finite()).reduce(T::max).unwrap();
    
    for y in y.iter_mut()
        .flatten()
    {
        *y = y.clamp(y_min, y_max)
    }
    
    let area = BitMapBackend::new(plot_path, PLOT_RES).into_drawing_area();
    
    area.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&area)
        .caption(plot_title, PLOT_CAPTION_FONT.into_font())
        .margin(PLOT_MARGIN)
        .x_label_area_size(PLOT_LABEL_AREA_SIZE)
        .y_label_area_size(PLOT_LABEL_AREA_SIZE)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    
    chart.configure_mesh()
        .set_all_tick_mark_size(0.1)
        .draw()?;
    
    for (i, (x, y)) in x.into_iter()
        .zip(y)
        .enumerate()
    {
        let color = Palette99::pick(i);
        chart.draw_series(LineSeries::new(
                x.into_iter()
                    .zip(y),
                &color
            ))?
        .label(format!("{}", i))
        .legend(move |(x, y)| Rectangle::new([(x + 5, y - 5), (x + 15, y + 5)], color.mix(0.5).filled()));
    }
    
    chart.configure_series_labels()
        .border_style(&BLACK)
        .draw()?;
        
    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present().expect("Unable to write result to file");

    Ok(())
}

pub fn plot_curve_2d<const NX: usize, const NY: usize>(
    plot_title: &str, plot_path: &str,
    x: [T; NX],
    y: [T; NY],
    f: impl Fn(T, T) -> T
) -> Result<(), Box<dyn std::error::Error>>
{
    use plotters::prelude::*;

    let area = SVGBackend::new(plot_path, PLOT_RES).into_drawing_area();
    
    let x_min = x.into_iter().reduce(T::min).unwrap();
    let x_max = x.into_iter().reduce(T::max).unwrap();
    
    let y_min = y.into_iter().reduce(T::min).unwrap();
    let y_max = y.into_iter().reduce(T::max).unwrap();

    let f_ref = &f;
    let f_values: Vec<T> = y.into_iter().flat_map(|y| x.into_iter().map(move |x| f_ref(x, y))).collect();

    let (z_min, z_max) = f_values.into_iter().map(|f| (f, f)).reduce(|a, b| (a.0.min(b.0), a.1.max(b.1))).unwrap();

    area.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&area)
        .caption(plot_title, PLOT_CAPTION_FONT)
        .set_all_label_area_size(PLOT_LABEL_AREA_SIZE)
        .build_cartesian_3d(x_min..x_max, z_min..z_max, y_min..y_max)?;

    chart.with_projection(isometric);
    
    chart.configure_axes()
        .light_grid_style(BLACK.mix(0.15))
        .max_light_lines(3)
        .draw()?;
    
    chart.draw_series(
            SurfaceSeries::xoz(
                x.into_iter(),
                y.into_iter(),
                f,
            )
            .style(BLUE.mix(0.2).filled()),
        )?
        .label("Surface")
        .legend(|(x, y)| Rectangle::new([(x + 5, y - 5), (x + 15, y + 5)], BLUE.mix(0.5).filled()));
    
    chart.configure_series_labels()
        .border_style(&BLACK)
        .draw()?;
    
    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present().expect("Unable to write result to file");

    Ok(())
}

pub fn plot_parametric_curve_2d<const NU: usize, const NV: usize>(
    plot_title: &str, plot_path: &str,
    u: [T; NU],
    v: [T; NV],
    f: impl Fn(T, T) -> [T; 3]
) -> Result<(), Box<dyn std::error::Error>>
{
    use plotters::prelude::*;

    let area = SVGBackend::new(plot_path, PLOT_RES).into_drawing_area();
    
    let f_ref = &f;
    let f_values: Vec<[T; 3]> = u.into_iter().flat_map(|u| v.into_iter().map(move |v| f_ref(u, v))).collect();

    let ([x_min, y_min, z_min], [x_max, y_max, z_max]) = f_values.into_iter()
        .map(|f| (f, f))
        .reduce(|a, b| (
            a.0.into_iter()
                .zip(b.0)
                .map(|(a, b)| a.min(b))
                .next_chunk()
                .unwrap(),
            a.1.into_iter()
                .zip(b.1)
                .map(|(a, b)| a.max(b))
                .next_chunk()
                .unwrap()
        )).unwrap();

    area.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&area)
        .caption(plot_title, PLOT_CAPTION_FONT)
        .set_all_label_area_size(PLOT_LABEL_AREA_SIZE)
        .build_cartesian_3d(x_min..x_max, z_min..z_max, y_min..y_max)?;

    chart.with_projection(isometric);
    
    chart.configure_axes()
        .light_grid_style(BLACK.mix(0.15))
        .max_light_lines(3)
        .draw()?;
    
    chart.draw_series(
            SurfaceSeries::xoz(
                u.into_iter(),
                v.into_iter(),
                f,
            )
            .map(|polygon| {
                let mut exception = false;
                let mut sum_theta = 0.0;
                let mut n_theta = 0;
                let points: Vec<(T, T, T)> = polygon.point_iter()
                    .into_iter()
                    .map(|&(u, [x, y, z], v)| {
                        let theta = v.atan2(u);
                        if theta.is_finite()
                        {
                            sum_theta += theta;
                            n_theta += 1;
                        }
                        if v == 0.0 && u < 0.0
                        {
                            exception = true;
                        }
                        (x, z, y)
                    })
                    .collect();
                let avg_theta = if exception {core::f32::consts::PI} else {sum_theta / n_theta as T};
                let c = (avg_theta as f64/core::f64::consts::TAU + 1.0) % 1.0;
                Polygon::new(points, HSLColor(c, 1.0, 0.5).mix(0.2).filled())
            })
        )?
        .label("Surface")
        .legend(|(x, y)| Rectangle::new([(x + 5, y - 5), (x + 15, y + 5)], BLUE.mix(0.5).filled()));
    
    chart.configure_series_labels()
        .border_style(&BLACK)
        .draw()?;
    
    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present().expect("Unable to write result to file");

    Ok(())
}

pub fn plot_curve_2d_rad<const NTHETA: usize, const NR: usize>(
    plot_title: &str, plot_path: &str,
    r: [T; NR],
    theta: [T; NTHETA],
    f: impl Fn(T, T) -> T
) -> Result<(), Box<dyn std::error::Error>>
where
    [(); 2*NR]:
{
    use plotters::prelude::*;

    let area = SVGBackend::new(plot_path, PLOT_RES).into_drawing_area();

    let r_max = r.into_iter().map(|r| r.abs()).reduce(T::max).unwrap();
    
    let theta_min = theta.into_iter().map(|theta| theta.abs()).reduce(T::min).unwrap();
    let theta_max = theta.into_iter().map(|theta| theta.abs()).reduce(T::max).unwrap();

    let f_ref = &f;
    let f_values: Vec<T> = r.into_iter().flat_map(|r| theta.into_iter().map(move |theta| f_ref(r, theta))).collect();
    let (z_min, z_max) = f_values.into_iter().map(|f| (f, f)).reduce(|a, b| (a.0.min(b.0), a.1.max(b.1))).unwrap();

    area.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&area)
        .caption(plot_title, PLOT_CAPTION_FONT)
        .set_all_label_area_size(PLOT_LABEL_AREA_SIZE)
        .build_cartesian_3d(-r_max..r_max, z_min..z_max, -r_max..r_max)?;

    chart.with_projection(isometric);
    
    chart.configure_axes()
        .light_grid_style(BLACK.mix(0.15))
        .max_light_lines(3)
        .draw()?;

    chart.draw_series(
            SurfaceSeries::xoz(
                r.into_iter(),
                theta.into_iter(),
                |r, theta| f(r, theta),
            )
            //.style_func(&|&c| HSLColor(c as f64, 1.0, 0.5).mix(0.2).filled())
            .map(|polygon| {
                let mut sum_theta = 0.0;
                let points: Vec<(T, T, T)> = polygon.point_iter()
                    .into_iter()
                    .map(|&(r, z, theta)| {sum_theta += theta; (r*theta.cos(), z, r*theta.sin())})
                    .collect();
                let avg_theta = sum_theta / points.len() as T;
                let c = (((avg_theta - theta_min)/(theta_max - theta_min)) as f64 + 1.0) % 1.0;
                Polygon::new(points, HSLColor(c, 1.0, 0.5).mix(0.2).filled())
            })
        )?
        .label("Radial surface")
        .legend(|(x, y)| Rectangle::new([(x + 5, y - 5), (x + 15, y + 5)], BLUE.mix(0.5).filled()));
    
    chart.configure_series_labels()
        .border_style(&BLACK)
        .draw()?;
    
    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present().expect("Unable to write result to file");

    Ok(())
}

pub fn plot_parametric_curve_2d_rad<const NU: usize, const NV: usize>(
    plot_title: &str, plot_path: &str,
    u: [T; NU],
    v: [T; NV],
    f: impl Fn(T, T) -> [T; 3]
) -> Result<(), Box<dyn std::error::Error>>
{
    use plotters::prelude::*;

    let area = SVGBackend::new(plot_path, PLOT_RES).into_drawing_area();

    let f_ref = &f;
    let f_values: Vec<[T; 3]> = u.into_iter().flat_map(|u| v.into_iter().map(move |v| f_ref(u, v))).collect();

    let ([_r_min, theta_min, z_min], [r_max, theta_max, z_max]) = f_values.into_iter()
        .map(|f| (f, f))
        .reduce(|a, b| (
            a.0.into_iter()
                .zip(b.0)
                .map(|(a, b)| a.min(b))
                .next_chunk()
                .unwrap(),
            a.1.into_iter()
                .zip(b.1)
                .map(|(a, b)| a.max(b))
                .next_chunk()
                .unwrap()
        ))
        .unwrap();

    area.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&area)
        .caption(plot_title, PLOT_CAPTION_FONT)
        .set_all_label_area_size(PLOT_LABEL_AREA_SIZE)
        .build_cartesian_3d(-r_max..r_max, z_min..z_max, -r_max..r_max)?;

    chart.with_projection(isometric);
    
    chart.configure_axes()
        .light_grid_style(BLACK.mix(0.15))
        .max_light_lines(3)
        .draw()?;

    chart.draw_series(
            SurfaceSeries::xoz(
                u.into_iter(),
                v.into_iter(),
                |u, v| f(u, v),
            )
            //.style_func(&|&c| HSLColor(c as f64, 1.0, 0.5).mix(0.2).filled())
            .map(|polygon| {
                let mut sum_theta = 0.0;
                let points: Vec<(T, T, T)> = polygon.point_iter()
                    .into_iter()
                    .map(|&(_, [r, theta, z], _)| {sum_theta += theta; (r*theta.cos(), z, r*theta.sin())})
                    .collect();
                let avg_theta = sum_theta / points.len() as T;
                let c = (((avg_theta - theta_min)/(theta_max - theta_min)) as f64 + 1.0) % 1.0;
                Polygon::new(points, HSLColor(c, 1.0, 0.5).mix(0.2).filled())
            })
        )?
        .label("Radial surface")
        .legend(|(x, y)| Rectangle::new([(x + 5, y - 5), (x + 15, y + 5)], BLUE.mix(0.5).filled()));
    
    chart.configure_series_labels()
        .border_style(&BLACK)
        .draw()?;
    
    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present().expect("Unable to write result to file");

    Ok(())
}